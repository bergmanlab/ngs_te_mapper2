#!/usr/bin/env python3

import sys
import os
import subprocess
import logging
import re
from datetime import datetime, timedelta
import pysam

"""
This script provides utility functions that are used by ngs_te_mapper program
"""


def format_time(time):
    d = datetime(1, 1, 1) + timedelta(seconds=time)
    if d.hour == 0 and d.minute == 0:
        return "%d seconds" % (d.second)
    elif d.hour == 0 and d.minute != 0:
        return "%d minutes %d seconds" % (d.minute, d.second)
    else:
        return "%d hours %d minutes %d seconds" % (d.hour, d.minute, d.second)


def create_soft_link(input, out_dir):
    link = os.path.join(out_dir, os.path.basename(input))
    if not os.path.isabs(input):
        input = os.path.abspath(input)
    if os.path.islink(link):
        os.remove(link)
    try:
        os.symlink(input, link)
    except Exception as e:
        print(e)
        logging.exception("Create symbolic link for " + input + " failed")
        sys.exit(1)
    return link


def parse_input(input_read, input_library, input_reference, out_dir):
    """
    Parse input files. If bam file is provided, convert to fasta format.
    """
    logging.info("Parsing input files...")

    # create symbolic link for the input file
    read = create_soft_link(input_read, out_dir)
    library = create_soft_link(input_library, out_dir)
    ref = create_soft_link(input_reference, out_dir)

    # unzip and merge input files, if muliple inputs were provided
    if ".gz" in input_read:
        fastq = out_dir + "/" + os.path.basename(input_read).replace(".gz", "")
        with open(fastq, "w") as output:
            subprocess.call(["gunzip", "-c", read], stdout=output)
    else:
        fastq = read

    return fastq, library, ref


def bam2fastq(bam, fastq):
    """
    Convert bam to fastq.
    """
    try:
        subprocess.call(["bedtools", "bamtofastq", "-i", bam, "-fq", fastq])
    except Exception as e:
        print(e)
        print("BAM to Fastq conversion failed, check input bam file, exiting...")
        sys.exit(1)


def get_bam_stats(bam, stat, ref_index, families, dir):
    # this script generates alignment stats
    bed_unique = dir + "/" + "unique.bed"
    bed_none = dir + "/" + "none.bed"
    with open(ref_index, "r") as input, open(bed_unique, "w") as unique, open(
        bed_none, "w"
    ) as none:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            out_line = "\t".join([entry[0], "1", entry[1]])
            if "chr" in line:
                unique.write(out_line + "\n")
            elif all(x not in line for x in families):
                none.write(out_line + "\n")

    bam_name = os.path.basename(bam).replace(".bam", "")
    sample_name = bam_name.split(".")[0]
    read_name = bam_name.split(".")[1]
    with open(stat, "a") as output:
        output.write(sample_name + "\t")
        output.write(read_name + "\t")
        num = subprocess.Popen(["samtools", "view", "-c", bam], stdout=subprocess.PIPE)
        total_reads = num.stdout.read().decode().replace("\n", "")
        output.write(total_reads + "\t")
        num = subprocess.Popen(
            ["samtools", "view", "-c", "-F", "260", bam], stdout=subprocess.PIPE
        )
        mapped_reads = num.stdout.read().decode().replace("\n", "")
        pt = "{:.1%}".format(int(mapped_reads) / int(total_reads))
        output.write(pt + "\t")
        # unique region
        readc, pt = count_reads(bam, bed_unique, total_reads)
        output.write("{:.1%}".format(pt) + "\t")
        pt_tes = 0
        # per family
        for family in families:
            bed = dir + "/" + family + ".bed"
            with open(ref_index, "r") as INDEX, open(bed, "w") as BED:
                for line in INDEX:
                    entry = line.replace("\n", "").split("\t")
                    out_line = "\t".join([entry[0], "1", entry[1]])
                    if family in line:
                        BED.write(out_line + "\n")
            readc, pt = count_reads(bam, bed, total_reads)
            pt_tes = pt_tes + pt
            output.write("{:.1%}".format(pt) + "\t")
            os.remove(bed)
        # focal TE regions
        output.write("{:.1%}".format(pt_tes) + "\t")
        # none region
        readc, pt = count_reads(bam, bed_none, total_reads)
        output.write("{:.1%}".format(pt) + "\t")
        output.write("\n")
        os.remove(bed_none)
        os.remove(bed_unique)


def get_lines(path):
    if os.path.isfile(path) == False or os.stat(path).st_size == 0:
        count = 0
    else:
        count = len(open(path).readlines())
    return count


def get_cluster(bam, bed, cutoff, window):
    # generate potential TE enriched clusters based on depth profile
    depth = bam + ".depth"
    with open(depth, "w") as output:
        subprocess.call(["samtools", "depth", bam, "-d", "0", "-Q", "1"], stdout=output)

    if os.path.isfile(depth) == False or os.stat(depth).st_size == 0:
        print("No depth info \n")
        return None

    depth_filter = depth + ".filter"
    with open(depth, "r") as input, open(depth_filter, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if int(entry[2]) >= cutoff:
                out_line = "\t".join(
                    [entry[0], str(entry[1]), str(int(entry[1]) + 1), str(entry[2])]
                )
                output.write(out_line + "\n")

    if os.path.isfile(depth_filter) == False or os.stat(depth_filter).st_size == 0:
        print("No depth info \n")
        return None

    bed_tmp = bed + ".tmp"
    with open(bed_tmp, "w") as output:
        subprocess.call(
            [
                "bedtools",
                "merge",
                "-d",
                str(window),
                "-c",
                "4",
                "-o",
                "mean",
                "-i",
                depth_filter,
            ],
            stdout=output,
        )

    # filter out non-chr entries
    with open(bed, "w") as output, open(bed_tmp, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if "chr" in entry[0]:
                output.write(line)

    os.remove(depth)
    os.remove(depth_filter)
    os.remove(bed_tmp)
    return None


def count_reads(bam, bed, total_reads):
    # cout number of reads mapped to a specific region
    command = "bedtools intersect -abam " + bam + " -b " + bed + " | samtools view -c"
    num = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    unique_reads = num.stdout.read().decode().replace("\n", "")
    pt = int(unique_reads) / int(total_reads)
    return unique_reads, pt


def file2set(file):
    file_set = set()
    with open(file, "r") as input:
        for line in input:
            file_set.add(line.replace("\n", ""))
    return file_set


def extract_reads(reads, read_file, out):
    """Extract reads from fasta using read ID list"""
    # read_ids = file2set(read_file)
    # record_dict = SeqIO.index(reads, "fastq")
    # with open(out, "wb") as output_handle:
    #     for key in read_ids:
    #         output_handle.write(record_dict.get_raw(key))

    # subset_fa = os.path.join(out, sample_name + ".subset.fa")

    command = "seqtk subseq " + reads + " " + read_file
    with open(out, "w") as output:
        subprocess.call(command, stdout=output, shell=True)


def mkdir(dir):
    try:
        os.mkdir(dir)
    except OSError:
        print("Creation of the directory %s failed" % dir)
    else:
        print("Successfully created the directory %s " % dir)


def sort_index_bam(bam, sorted_bam, thread):
    """
    Sort and index bam file
    """
    try:
        subprocess.call(["samtools", "sort", "-@", str(thread), "-o", sorted_bam, bam])
        subprocess.call(["samtools", "index", "-@", str(thread), sorted_bam])
    except Exception as e:
        print(e)
        print("Sort and index BAM file failed, exiting...")
        sys.exit(1)


def make_bam(fq, ref, thread, bam, mapper="minimap2"):
    # alignment and generate sorted bam file
    sam = bam + ".sam"

    if mapper == "minimap2":
        presets = "sr"
        with open(sam, "w") as output:
            subprocess.call(
                ["minimap2", "-ax", presets, "-v", "0", "-t", str(thread), ref, fq],
                stdout=output,
            )
    else:
        with open(sam, "w") as output:
            subprocess.call(
                ["bwa", "mem", "-v", "0", "-t", str(thread), ref, fq], stdout=output
            )

    command = (
        "samtools view -Sb -t "
        + ref
        + " "
        + sam
        + " | "
        + "samtools sort -@ "
        + str(thread)
        + " -o "
        + bam
    )
    subprocess.call(command, shell=True)
    subprocess.call(["samtools", "index", bam])
    os.remove(sam)


def subset_bam_ids(bam_in, bam_out, contigs, ids, thread):
    # subset bam file using read IDs
    # # https://bioinformatics.stackexchange.com/questions/3380/how-to-subset-a-bam-by-a-list-of-qnames
    bam_tmp = bam_out + ".tmp"
    with open(bam_tmp, "w") as output:
        command = (
            "samtools view -h "
            + bam_in
            + " "
            + contigs
            + " | awk 'FNR==NR {reads[$1]||$0 ~ /@/;next} /^@/||($1 in reads) {print $0}' "
            + ids
            + " - "
            # + " | awk '{if($0 ~ /^@/ || $6 ~ /S/ || $6 ~ /H/) {print $0}}' "
            + " | samtools view -Sb -"
        )
        subprocess.call(command, shell=True, stdout=output)
    sort_index_bam(bam_tmp, bam_out, thread)


def repeatmask(ref, library, outdir, thread, augment=False):
    try:
        subprocess.call(
            [
                "RepeatMasker",
                "-dir",
                outdir,
                "-gff",
                "-s",
                "-nolow",
                "-no_is",
                "-xsmall",
                "-e",
                "ncbi",
                "-lib",
                library,
                "-pa",
                str(thread),
                ref,
            ]
        )
        ref_rm = os.path.join(outdir, os.path.basename(ref) + ".masked")
        open(ref_rm, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking failed, exiting...")
        sys.exit(1)
    if augment:
        ref_rm_aug = ref_rm + ".aug"
        subprocess.call(["cat", ref_rm, library], stdout=output)
        return ref_rm_aug
    else:
        return ref_rm


def get_family_bam(bam_in, bam_out, family, thread):
    # filter bam file for reads partially mapped to a given TE family
    bam_tmp = bam_out + ".tmp"
    with open(bam_tmp, "w") as output:
        command = (
            "samtools view -h "
            + bam_in
            + " "
            + family
            + " | awk '{if($0 ~ /^@/ || $6 ~ /S/ || $6 ~ /H/) {print $0}}' "
            + " | samtools view -Sb -"
        )
        subprocess.call(command, shell=True, stdout=output)

    # sort and index
    sort_index_bam(bam_tmp, bam_out, thread)
    os.remove(bam_tmp)
    # pysam.sort("-o", bam_out, bam_tmp)
    # pysam.index(bam_out)


def get_bam_id(bam_in, pattern, outfile):
    bam_open = pysam.AlignmentFile(bam_in, mode="r")

    with open(outfile, "w") as output:
        for read in bam_open.fetch():
            read = str(read)
            read_id = read.split("\t")[0]
            cigar = read.split("\t")[5].replace("H", "S")
            cigar_pattern = "".join(re.findall("[SM]+", cigar))
            if cigar_pattern == pattern:
                output.write(read_id + "\n")


def filter_bam(bam_in, read_list, bam_out):
    # filter bam file based on provided read list
    command = (
        "picard FilterSamReads I="
        + bam_in
        + " O="
        + bam_out
        + " READ_LIST_FILE="
        + read_list
        + " FILTER=includeReadList"
    )
    with open(bam_out, "w") as output:
        subprocess.call(command, stdout=output, shell=True)
    subprocess.call(["samtools", "index", bam_out])


def get_family_bed(args):
    family = args[0]
    bam = args[1]
    ref_rm = args[2]
    outdir = args[3]
    mapper = args[4]
    contigs = args[5]
    ngs_te_mapper = args[6]

    family_dir = os.path.join(outdir, family)
    mkdir(family_dir)
    # get CIGAR alignments to TE family
    family_bam = family_dir + "/" + family + ".te.cigar.bam"
    get_family_bam(bam_in=bam, bam_out=family_bam, family=family, thread=1)

    sm_ids = family_bam.replace(".bam", ".sm.ids")
    get_bam_id(bam_in=family_bam, pattern="SM", outfile=sm_ids)
    ms_ids = family_bam.replace(".bam", ".ms.ids")
    get_bam_id(bam_in=family_bam, pattern="MS", outfile=ms_ids)

    # step three: reads alignment to reference genome
    ms_bam = family_dir + "/" + family + ".ref.cigar.ms.bam"
    sm_bam = family_dir + "/" + family + ".ref.cigar.sm.bam"
    if ngs_te_mapper:
        cigar_fq = family_dir + "/" + family + ".cigar.fastq"
        bam2fastq(family_bam, cigar_fq)
        # need to realign to TE cigar reads to masked ref
        sm_fq = family_dir + "/" + family + ".te.cigar.sm.fastq"
        extract_reads(cigar_fq, sm_ids, sm_fq)
        ms_fq = family_dir + "/" + family + ".te.cigar.ms.fastq"
        extract_reads(cigar_fq, ms_ids, ms_fq)

        # map reads to masked reference genome
        make_bam(fq=sm_fq, ref=ref_rm, thread=1, bam=sm_bam, mapper=mapper)
        make_bam(fq=ms_fq, ref=ref_rm, thread=1, bam=ms_bam, mapper=mapper)
    else:
        # subset BAM using read IDs using awk way
        subset_bam_ids(
            bam_in=bam, bam_out=ms_bam, contigs=contigs, ids=ms_ids, thread=1
        )
        subset_bam_ids(
            bam_in=bam, bam_out=sm_bam, contigs=contigs, ids=sm_ids, thread=1
        )

    # get insertion candidate using coverage profile
    sm_bed = family_dir + "/" + family + ".sm.bed"
    get_cluster(sm_bam, sm_bed, cutoff=1, window=100)

    ms_bed = family_dir + "/" + family + ".ms.bed"
    get_cluster(ms_bam, ms_bed, cutoff=1, window=100)

    # get insertion candidate
    get_nonref(sm_bed, ms_bed, family_dir, family)


def get_nonref(bed1, bed2, outdir, family, tsd_max=20):
    overlap = outdir + "/" + family + ".overlap.tsv"
    if os.path.isfile(bed1) and os.path.isfile(bed2):
        with open(overlap, "w") as output:
            subprocess.call(
                ["bedtools", "intersect", "-wao", "-a", bed1, "-b", bed2], stdout=output
            )

    # parse overlap
    if os.path.isfile(overlap):
        nonref = outdir + "/" + family + ".nonref.bed"
        with open(overlap, "r") as input, open(nonref, "w") as output:
            for line in input:
                entry = line.replace("\n", "").split("\t")
                if entry[4] != "." and int(entry[8]) <= tsd_max:
                    chr = entry[0]
                    if abs(int(entry[1]) - int(entry[6])) < abs(
                        int(entry[2]) - int(entry[5])
                    ):
                        start = entry[1]
                        end = entry[6]
                    else:
                        start = entry[5]
                        end = entry[2]
                    score = (float(entry[3]) + float(entry[7])) / 2
                    score = "{:.2f}".format(score)
                    strand = "."
                    out_line = "\t".join(
                        [chr, str(start), str(end), family, str(score), strand]
                    )
                    output.write(out_line + "\n")
        # os.remove(overlap)

        # if os.path.isfile(nonref):
        #     nonref_norm = sample_dir + "/" + sample_name + "." + family + ".nonref.norm.bed"
        #     with open(nonref_norm, "w") as output:
        #         subprocess.call(["bedtools", "intersect", "-a", nonref, "-b", norm, "-u"], stdout=output)
        #     nonref_norm_count = get_lines(nonref_norm)
        # else:
        #     nonref_norm_count = 0
    # else:
    #     nonref_count = 0
    #     nonref_norm_count = 0

    # with open(stats, "a") as output:
    #     output.write(sample_name+"\t")
    #     output.write(family+"\t")
    #     output.write("non-ref\t")
    #     output.write(str(nonref_count)+"\t")
    #     output.write(str(nonref_norm_count)+"\n")


def merge_bed(bed_in, bed_out):
    # merge bed files from all families, check overlap within and betweeen families, merge or remove entried if necessary
    bed_out_tmp = bed_out + ".tmp"
    with open(bed_out_tmp, "w") as output:
        for bed in bed_in:
            if os.path.isfile(bed) and os.stat(bed).st_size != 0:
                with open(bed, "r") as input:
                    for line in input:
                        output.write(line)

    # sort bed files
    with open(bed_out, "w") as output:
        subprocess.call(["bedtools", "sort", "-i", bed_out_tmp], stdout=output)
    os.remove(bed_out_tmp)


def get_ref_te(
    bed1, bed2, norm, gff, family, stat, sample_dir, tmp_dir, sample_name, window=100
):
    # get reference count
    ref_te = tmp_dir + "/" + family + ".ref.bed"
    with open(ref_te, "w") as output, open(gff, "r") as input:
        for line in input:
            check_family = "Name=" + family + "{}"
            if check_family in line:
                entry = line.replace("\n", "").split("\t")
                family = entry[8].split(";")[1]
                family = re.sub("Name=", "", family)
                family = re.sub("{}.*", "", family)
                out_line = "\t".join(
                    [entry[0], entry[3], entry[4], family, ".", entry[6]]
                )
                output.write(out_line + "\n")
    ref_count = get_lines(ref_te)
    ref_te_norm = tmp_dir + "/" + family + ".ref.norm.bed"
    with open(ref_te_norm, "w") as output:
        subprocess.call(
            ["bedtools", "intersect", "-a", ref_te, "-b", norm, "-u"], stdout=output
        )

    # calculate clusters that jointly support ref TEs (all, norm) with a percentage
    ref_set1 = bed1.replace(".bed", ".ref.bed")
    ref_set1_norm = bed1.replace(".bed", ".ref.norm.bed")
    if os.path.isfile(bed1):
        with open(ref_set1, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te,
                    "-b",
                    bed1,
                    "-u",
                ],
                stdout=output,
            )
        ref_set1_count = get_lines(ref_set1)

        with open(ref_set1_norm, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te_norm,
                    "-b",
                    bed1,
                    "-u",
                ],
                stdout=output,
            )
        ref_set1_norm_count = get_lines(ref_set1_norm)
    else:
        ref_set1_count = 0
        ref_set1_norm_count = 0

    ref_set2 = bed2.replace(".bed", ".ref.bed")
    ref_set2_norm = bed2.replace(".bed", ".ref.norm.bed")
    if os.path.isfile(bed2):
        with open(ref_set2, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te,
                    "-b",
                    bed2,
                    "-u",
                ],
                stdout=output,
            )
        ref_set2_count = get_lines(ref_set2)

        with open(ref_set2_norm, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    ref_te_norm,
                    "-b",
                    bed2,
                    "-u",
                ],
                stdout=output,
            )
        ref_set2_norm_count = get_lines(ref_set2_norm)
    else:
        ref_set2_count = 0
        ref_set2_norm_count = 0

    # joinly support ref
    if os.path.isfile(ref_set1) and os.path.isfile(ref_set2):
        ref_both = tmp_dir + "/" + sample_name + "." + family + ".ref.bed"
        with open(ref_both, "w") as output:
            subprocess.call(
                ["bedtools", "intersect", "-a", ref_set1, "-b", ref_set2, "-u"],
                stdout=output,
            )
        ref_both_count = get_lines(ref_both)
    else:
        ref_both_count = 0

    if os.path.isfile(ref_set1_norm) and os.path.isfile(ref_set2_norm):
        ref_both_norm = sample_dir + "/" + sample_name + "." + family + ".ref.norm.bed"
        with open(ref_both_norm, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "intersect",
                    "-a",
                    ref_set1_norm,
                    "-b",
                    ref_set2_norm,
                    "-u",
                ],
                stdout=output,
            )
        ref_both_norm_count = get_lines(ref_both_norm)
    else:
        ref_both_norm_count = 0

    # write stats to summary
    with open(stat, "a") as output:
        output.write(sample_name + "\t")
        output.write(family + "\t")
        output.write("ref\t")
        out_num = (
            str(ref_both_count)
            + "("
            + str(ref_set1_count)
            + ","
            + str(ref_set2_count)
            + ")"
        )
        output.write(out_num + "\t")
        out_num = (
            str(ref_both_norm_count)
            + "("
            + str(ref_set1_norm_count)
            + ","
            + str(ref_set2_norm_count)
            + ")"
        )
        output.write(out_num + "\n")

    # clean tmp files
    os.remove(ref_te)
    os.remove(ref_te_norm)