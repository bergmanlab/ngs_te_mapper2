#!/usr/bin/env python3

import sys
import os
import subprocess
import logging
import re
from datetime import datetime, timedelta
from statistics import mean
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


def parse_input(input_reads, input_library, input_reference, out_dir):
    """
    Parse input files. If bam file is provided, convert to fasta format.
    """
    logging.info("Parsing input files...")

    # create symbolic link for the input file
    library = create_soft_link(input_library, out_dir)
    ref = create_soft_link(input_reference, out_dir)

    input_reads = input_reads.replace(" ", "").split(",")
    # unzip and merge input files, if muliple inputs were provided
    if len(input_reads) == 0:
        logging.exception("no reads are provided, check your input files")
        sys.exit(1)
    reads_copy = []
    for read in input_reads:
        reads_copy.append(create_soft_link(read, out_dir))
    if len(reads_copy) == 1:
        read = reads_copy[0]
        prefix = get_prefix(read)
        if ".gz" in read:
            fastq = read.replace(".gz", "")
            with open(fastq, "w") as output:
                subprocess.call(["gunzip", "-c", read], stdout=output)
        else:
            fastq = read
    else:
        prefix = get_prefix(reads_copy[0])
        fastq = os.path.join(out_dir, prefix + ".fastq")
        with open(fastq, "w") as output:
            for read in reads_copy:
                if ".gz" in read:
                    subprocess.call(["gunzip", "-c", read], stdout=output)
                else:
                    subprocess.call(["cat", read], stdout=output)
    return prefix, fastq, library, ref


def get_prefix(path):
    no_path = os.path.basename(path)
    if no_path[-3:] == ".gz":
        no_path = no_path[:-3]
    no_ext = ".".join(no_path.split(".")[:-1])
    return no_ext


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
                    [entry[0], str(int(entry[1]) - 1), str(entry[1]), str(entry[2])]
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

    # filter out non-chr entries convert to half close coordinate (bed format)
    with open(bed, "w") as output, open(bed_tmp, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if "chr" in entry[0]:
                out_line = "\t".join(
                    [entry[0], entry[1], str(int(entry[2]) + 1), entry[3]]
                )
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


def make_bam(fq, ref, thread, bam, mapper="bwa"):
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
        gff = os.path.join(outdir, os.path.basename(ref) + ".out.gff")
        if not os.path.isfile(ref_rm):
            ref_rm_out = os.path.join(outdir, os.path.basename(ref) + ".out")
            with open(ref_rm_out, "r") as input:
                for line in input:
                    if "There were no repetitive sequences detected" in line:
                        print("No repetitive sequences detected")
                        ref_rm = ref
                        gff = None
                        rm_bed = None
                    else:
                        raise Exception("Repeatmasking failed, exiting...")
        else:
            rm_bed = os.path.join(outdir, os.path.basename(ref) + ".bed")
            parse_rm_out(gff, rm_bed)
            open(ref_rm, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking failed, exiting...")
        sys.exit(1)
    if augment:
        ref_rm_aug = ref_rm + ".aug"
        with open(ref_rm_aug, "w") as output:
            subprocess.call(["cat", ref_rm, library], stdout=output)
        return ref_rm_aug, rm_bed
    else:
        return ref_rm, rm_bed


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
    rm_bed = args[3]
    outdir = args[4]
    mapper = args[5]
    contigs = args[6]
    experiment = args[7]
    tsd_max = args[8]
    gap_max = args[9]
    window = args[10]

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
    if not experiment:
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
    get_cluster(
        sm_bam, sm_bed, cutoff=1, window=window
    )  # TODO: add this option to args?

    sm_bed_refined = family_dir + "/" + family + ".refined.sm.bed"
    if os.path.isfile(sm_bed) and os.stat(sm_bed).st_size != 0:
        refine_breakpoint(sm_bed_refined, sm_bed, sm_bam, "SM")

    ms_bed = family_dir + "/" + family + ".ms.bed"
    get_cluster(ms_bam, ms_bed, cutoff=1, window=window)

    ms_bed_refined = family_dir + "/" + family + ".refined.ms.bed"
    if os.path.isfile(ms_bed) and os.stat(ms_bed).st_size != 0:
        refine_breakpoint(ms_bed_refined, ms_bed, ms_bam, "MS")

    if (
        rm_bed is not None
        and os.path.isfile(sm_bed_refined)
        and os.path.isfile(ms_bed_refined)
    ):
        get_ref(sm_bed_refined, ms_bed_refined, rm_bed, family_dir, family)
        # sm_bed_filtered = family_dir + "/" + family + ".filtered.sm.bed"
        # ref_te_dict, sm_bed_filtered = get_ref(
        #     sm_bed_refined, rm_bed, family_dir, family
        # )

    # get insertion candidate
    if os.path.isfile(sm_bed_refined) and os.path.isfile(ms_bed_refined):
        get_nonref(sm_bed_refined, ms_bed_refined, family_dir, family, tsd_max, gap_max)


def get_genome_file(ref):
    subprocess.call(["samtools", "faidx", ref])
    ref_index = ref + ".fai"
    genome = ref + ".genome"
    with open(ref_index, "r") as input, open(genome, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            out_line = "\t".join([entry[0], entry[1]])
            output.write(out_line + "\n")
    return genome


def get_af(
    bed, ref, reads, thread, out_dir, sample_prefix, method="max", slop=150, offset=10
):
    # create masked genome (only allow bases around non-ref TEs)
    genome = get_genome_file(ref)
    expand_bed = out_dir + "/" + sample_prefix + ".nonref.expand.bed"
    with open(expand_bed, "w") as output:
        subprocess.call(
            ["bedtools", "slop", "-i", bed, "-g", genome, "-b", str(slop)],
            stdout=output,
        )

    complement_bed = out_dir + "/" + sample_prefix + ".nonref.complement.bed"
    with open(complement_bed, "w") as output:
        subprocess.call(
            ["bedtools", "complement", "-i", expand_bed, "-g", genome], stdout=output
        )

    masked_ref = ref + ".nonref.masked"
    subprocess.call(
        ["bedtools", "maskfasta", "-fi", ref, "-bed", complement_bed, "-fo", masked_ref]
    )
    subprocess.call(["bwa", "index", masked_ref])
    # map raw reads to masked genome
    bam = out_dir + "/" + sample_prefix + ".nonref.bam"
    make_bam(fq=reads, ref=masked_ref, thread=thread, bam=bam)

    # for each site, figure out AF
    af_bed = bed.replace(".bed", ".af.bed")
    samfile = pysam.AlignmentFile(bam, "rb")
    with open(bed, "r") as input, open(af_bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chromosome = entry[0]
            start = int(entry[1])
            end = int(entry[2])
            info = entry[3]
            score = entry[4]
            strand = entry[5]

            n_cigar1 = count_read(samfile, chromosome, start - offset, start, "MS")
            n_cigar2 = count_read(samfile, chromosome, end, end + offset, "SM")

            n_ref = count_read(samfile, chromosome, start, end, "M")

            af1 = n_cigar1 / (n_cigar1 + n_ref)
            af2 = n_cigar2 / (n_cigar2 + n_ref)

            if method == "max":
                af = round(max(af1, af2), 2)
            else:
                af = round(mean([af1, af2]), 2)

            # output
            info_new = "|".join(
                [info, str(af), str(n_cigar1), str(n_cigar2), str(n_ref)]
            )
            out_line = "\t".join(
                [chromosome, str(start), str(end), info_new, score, strand]
            )
            output.write(out_line + "\n")
    return af_bed


def count_read(samfile, chromosome, start, end, cigar_pattern):
    n_read = 0
    for read in samfile.fetch(chromosome, start, end):
        pattern, offset = parse_cigar(read.cigar)
        if pattern == cigar_pattern:
            n_read = n_read + 1
    return n_read


def parse_cigar(cigar):
    offset = 0
    pattern = ""
    for item in cigar:
        if item[0] == 0:
            offset = offset + item[1]
            if "M" not in pattern:
                pattern = pattern + "M"
        elif item[0] == 1:
            continue
        elif item[0] == 2:
            offset = offset + item[1]
        elif item[0] == 4 or item[0] == 5:
            pattern = pattern + "S"
        else:
            print("I don't recognize this string:" + str(item[0]) + "\n")
    return pattern, offset


def refine_breakpoint(new_bed, old_bed, bam, cigar_type):
    samfile = pysam.AlignmentFile(bam, "rb")
    with open(old_bed, "r") as input, open(new_bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chromosome = entry[0]
            cov_start = int(entry[1])
            cov_end = int(entry[2])
            cov_count = entry[3]
            breakpoints = dict()
            cigars = {"SM": 0, "MS": 0}
            for read in samfile.fetch(chromosome, cov_start, cov_end):
                ref_start = read.reference_start
                pattern, offset = parse_cigar(read.cigar)
                # if offset > 150 or offset < -150:
                #     print("strange offset:" + str(offset) + "\n")
                #     print(read.cigar)
                #     print(read.reference_start)
                #     print(read.query_name)
                if pattern == "SM":
                    cigars["SM"] = cigars["SM"] + 1
                    breakpoint = ref_start
                elif pattern == "MS":
                    cigars["MS"] = cigars["MS"] + 1
                    breakpoint = ref_start + offset
                else:
                    # print("weird pattern:" + pattern)
                    # print(read.query_name)
                    # print(read.cigar)
                    breakpoint = None
                if breakpoint:
                    if breakpoint in breakpoints:
                        breakpoints[breakpoint] = breakpoints[breakpoint] + 1
                    else:
                        breakpoints[breakpoint] = 1

            cigar_joint = get_keys(cigars)
            breakpoint_joint = get_keys(breakpoints)

            if cigar_joint:
                if len(cigar_joint) > 1:
                    cigar_joint = "mix"
                else:
                    cigar_joint = cigar_joint[0]
            if breakpoint_joint:
                breakpoint_joint = round(mean(breakpoint_joint))

            if cigar_joint == "MS":
                refined_start = cov_start
                refined_end = breakpoint_joint
            elif cigar_joint == "SM":
                refined_start = breakpoint_joint
                refined_end = cov_end
            else:
                refined_start = cov_start
                refined_end = cov_end

            if (cigar_joint == "SM" or cigar_joint == "MS") and breakpoint_joint:
                cigar_count = cigars[cigar_joint]
                out_line = "\t".join(
                    [
                        chromosome,
                        str(refined_start),
                        str(refined_end),
                        str(cigar_count),
                    ]
                )
                output.write(out_line + "\n")
    # sys.exit(1)


def get_keys(dictionary):
    if bool(dictionary):
        max_value = max(dictionary.values())
        indices = [i for i, x in dictionary.items() if x == max_value]
        return indices
    else:
        return None


def get_nonref(bed1, bed2, outdir, family, tsd_max, gap_max):
    overlap = outdir + "/" + family + ".overlap.tsv"
    if os.path.isfile(bed1) and os.path.isfile(bed2):
        with open(overlap, "w") as output:
            # subprocess.call(
            #     ["bedtools", "intersect", "-wao", "-a", bed1, "-b", bed2], stdout=output
            # )
            subprocess.call(
                ["bedtools", "window", "-w", str(gap_max), "-a", bed1, "-b", bed2],
                stdout=output,
            )

    # parse overlap
    if os.path.isfile(overlap) and os.stat(overlap).st_size != 0:
        nonref = outdir + "/" + family + ".nonref.bed"
        with open(overlap, "r") as input, open(nonref, "w") as output:
            for line in input:
                ins_pass = True
                entry = line.replace("\n", "").split("\t")
                chromosome = entry[0]
                if (
                    (int(entry[1]) - int(entry[5])) > 0
                    and (int(entry[2]) - int(entry[6])) > 0
                ) or (
                    (int(entry[5]) - int(entry[1])) > 0
                    and (int(entry[6]) - int(entry[2])) > 0
                ):  # get rid of entries if one is within another
                    if abs(int(entry[1]) - int(entry[6])) < abs(
                        int(entry[2]) - int(entry[5])
                    ):
                        if int(entry[1]) < int(entry[6]):
                            start = entry[1]
                            end = entry[6]
                        else:
                            start = entry[6]
                            end = entry[1]
                        strand = "-"
                        dist = int(entry[1]) - int(
                            entry[6]
                        )  # positive: gap; negative: overlap
                    else:
                        if int(entry[2]) < int(entry[5]):
                            start = entry[2]
                            end = entry[5]
                        else:
                            start = entry[5]
                            end = entry[2]
                        strand = "+"
                        dist = int(entry[5]) - int(
                            entry[2]
                        )  # positive: gap; negative: overlap
                    if dist < 0:
                        if -dist > tsd_max:
                            ins_pass = False
                    else:
                        if dist > gap_max:
                            ins_pass = False

                    if ins_pass:
                        # score = (float(entry[3]) + float(entry[7])) / 2
                        # score = "{:.2f}".format(score)
                        # score = "|".join([entry[3], entry[7]])
                        score = "."
                        family_info = "|".join([family, str(dist)])
                        out_line = "\t".join(
                            [
                                chromosome,
                                str(start),
                                str(end),
                                family_info,
                                str(score),
                                strand,
                            ]
                        )
                        output.write(out_line + "\n")
        # os.remove(overlap)


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


def parse_rm_out(rm_gff, bed):
    with open(bed, "w") as output, open(rm_gff, "r") as input:
        for line in input:
            if "RepeatMasker" in line:
                entry = line.replace("\n", "").split("\t")
                family = entry[8].split(" ")[1]
                family = re.sub('"Motif:', "", family)
                family = re.sub('"', "", family)
                out_line = "\t".join(
                    [entry[0], entry[3], entry[4], family, ".", entry[6]]
                )
                output.write(out_line + "\n")


def get_ref(cov_bed, rm_bed, out_dir, family, window=50):
    # calculate clusters that jointly support ref TEs (all, norm) with a percentage
    ref_te_cov = cov_bed.replace(".bed", ".ref.cov.bed")
    if os.path.isfile(cov_bed):
        with open(ref_te_cov, "w") as output:
            subprocess.call(
                [
                    "bedtools",
                    "window",
                    "-w",
                    str(window),
                    "-a",
                    rm_bed,
                    "-b",
                    cov_bed,
                    "-u",
                ],
                stdout=output,
            )
    # with open(ref_te_cov, "r") as input:
    #     for line in input:

    # ref_ms = bed2.replace(".bed", ".ref.bed")
    # if os.path.isfile(bed2):
    #     with open(ref_ms, "w") as output:
    #         subprocess.call(
    #             [
    #                 "bedtools",
    #                 "window",
    #                 "-w",
    #                 str(window),
    #                 "-a",
    #                 rm_bed,
    #                 "-b",
    #                 bed2,
    #                 "-u",
    #             ],
    #             stdout=output,
    #         )
    # go over those output files and build dict for reference TEs

    # if os.path.isfile(ref_sm) and os.path.isfile(ref_ms):
    #     ref_both = out_dir + "/" + family + ".ref.bed"
    #     with open(ref_both, "w") as output:
    #         subprocess.call(
    #             ["bedtools", "intersect", "-a", ref_sm, "-b", ref_ms, "-u"],
    #             stdout=output,
    #         )