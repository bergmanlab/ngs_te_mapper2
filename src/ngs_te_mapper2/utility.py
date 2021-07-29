#!/usr/bin/env python3

import sys
import os
import subprocess
import logging
import re
from Bio import SeqIO
import math
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


def fix_fasta_header(fasta_input, fasta_output):
    sequences = []
    with open(fasta_input, "r") as input, open(fasta_output, "w") as output:
        for record in SeqIO.parse(input, "fasta"):
            seq_name = str(record.id)
            if "#" in seq_name:
                seq_name = seq_name.split("#")[0]
                record.id = seq_name
                record.description = ""
            sequences.append(record)
        SeqIO.write(sequences, output, "fasta")


def parse_input(input_reads, input_library, input_reference, out_dir):
    """
    Parse input files.
    """
    logging.info("Parsing input files...")

    # create symbolic link for the input file
    library = create_soft_link(input_library, out_dir)
    library_fix = library + ".fix"
    fix_fasta_header(library, library_fix)

    ref = create_soft_link(input_reference, out_dir)
    ref_fix = ref + ".fix"
    fix_fasta_header(ref, ref_fix)

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
        prefix = get_prefix(read).replace("+", "plus")
        if ".gz" in read:
            fastq = read.replace(".gz", "")
            with open(fastq, "w") as output:
                subprocess.call(["gunzip", "-c", read], stdout=output)
        else:
            fastq = read
    else:
        prefix = get_prefix(reads_copy[0]).replace("+", "plus")
        fastq = os.path.join(out_dir, prefix + ".merge.fastq")
        with open(fastq, "w") as output:
            for read in reads_copy:
                if ".gz" in read:
                    subprocess.call(["gunzip", "-c", read], stdout=output)
                else:
                    subprocess.call(["cat", read], stdout=output)
    return prefix, fastq, library_fix, ref_fix


def get_masked_genome(ref, outdir, bed):
    ref_masked = os.path.join(outdir, os.path.basename(ref) + ".masked")
    subprocess.call(
        ["bedtools", "maskfasta", "-fi", ref, "-bed", bed, "-fo", ref_masked]
    )
    return ref_masked


def get_prefix(path):
    """
    Get prefix of the input fastq files.
    """
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


def get_lines(path):
    if os.path.isfile(path) == False or os.stat(path).st_size == 0:
        count = 0
    else:
        count = len(open(path).readlines())
    return count


def get_cluster(bam, bed, cutoff, window, family):
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
            contig = entry[0]
            if contig != family:
                # if "chr" in entry[0]:
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


def repeatmask(ref, library, outdir, thread):
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
        gff3 = os.path.join(outdir, os.path.basename(ref) + ".out.gff3")
        if not os.path.isfile(ref_rm):
            ref_rm_out = os.path.join(outdir, os.path.basename(ref) + ".out")
            with open(ref_rm_out, "r") as input:
                for line in input:
                    if "There were no repetitive sequences detected" in line:
                        print("No repetitive sequences detected")
                        ref_rm = ref
                        gff = None
                        gff3 = None
                    else:
                        raise Exception("Repeatmasking failed, exiting...")
        else:
            parse_rm_out(gff, gff3)
            open(ref_rm, "r")
    except Exception as e:
        print(e)
        print("Repeatmasking failed, exiting...")
        sys.exit(1)
    return ref_rm, gff3


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
    tsd_max = args[7]
    gap_max = args[8]
    window = args[9]
    min_mapq = args[10]

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

    # get insertion candidate using coverage profile
    sm_bed = family_dir + "/" + family + ".sm.bed"
    get_cluster(
        sm_bam, sm_bed, cutoff=1, window=window, family=family
    )  # TODO: add this option to args?

    sm_bed_refined = family_dir + "/" + family + ".refined.sm.bed"
    if os.path.isfile(sm_bed) and os.stat(sm_bed).st_size != 0:
        refine_breakpoint(sm_bed_refined, sm_bed, sm_bam, "SM", min_mapq)

    sm_bed_unique = family_dir + "/" + family + ".unique.sm.bed"
    if os.path.isfile(sm_bed_refined) and os.stat(sm_bed_refined).st_size != 0:
        bed_rm_dup(sm_bed_refined, sm_bed_unique)

    ms_bed = family_dir + "/" + family + ".ms.bed"
    get_cluster(ms_bam, ms_bed, cutoff=1, window=window, family=family)

    ms_bed_refined = family_dir + "/" + family + ".refined.ms.bed"
    if os.path.isfile(ms_bed) and os.stat(ms_bed).st_size != 0:
        refine_breakpoint(ms_bed_refined, ms_bed, ms_bam, "MS", min_mapq)

    ms_bed_unique = family_dir + "/" + family + ".unique.ms.bed"
    if os.path.isfile(ms_bed_refined) and os.stat(ms_bed_refined).st_size != 0:
        bed_rm_dup(ms_bed_refined, ms_bed_unique)

    if (
        rm_bed is not None
        and os.path.isfile(sm_bed_unique)
        and os.path.isfile(ms_bed_unique)
    ):
        get_ref(sm_bed_unique, ms_bed_unique, rm_bed, family_dir, family)

    # get insertion candidate
    if os.path.isfile(sm_bed_unique) and os.path.isfile(ms_bed_unique):
        get_nonref(sm_bed_unique, ms_bed_unique, family_dir, family, tsd_max, gap_max)


def bed_rm_overlap(bed_in, bed_out):
    bed_merge = bed_in + ".redundant"
    with open(bed_merge, "w") as output:
        command = 'bedtools merge -d 0 -o collapse -c 2,3,4,5,6 -delim "," -i ' + bed_in
        subprocess.call(command, shell=True, stdout=output)

    with open(bed_merge, "r") as input, open(bed_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if "," not in entry[3]:
                chromosome = entry[0]
                start = entry[3]
                end = entry[4]
                info = entry[5]
                score = entry[6]
                strand = entry[7]
                out_line = "\t".join([chromosome, start, end, info, score, strand])
                output.write(out_line + "\n")
    os.remove(bed_merge)


def bed_rm_duplicate(bed_in, bed_out):
    with open(bed_out, "w") as output:
        command = "cat " + bed_in + " | sort | uniq"
        subprocess.call(command, shell=True, stdout=output)


def bed_rm_dup(bed_in, bed_out):
    bed_merge = bed_in + ".merge.tmp"
    with open(bed_merge, "w") as output:
        command = 'bedtools merge -d 0 -o collapse -c 2,3,4,5,6 -delim "," -i ' + bed_in
        subprocess.call(command, shell=True, stdout=output)

    with open(bed_merge, "r") as input, open(bed_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chromosome = entry[0]
            if "," in entry[3]:
                score = entry[5].split(",")
                score = [int(i) for i in score]
                idx = score.index(max(score))
                start = entry[3].split(",")[idx]
                end = entry[4].split(",")[idx]
                score = entry[5].split(",")[idx]
                cigar = entry[6].split(",")[idx]
                bp = entry[7].split(",")[idx]
            else:
                start = entry[3]
                end = entry[4]
                score = entry[5]
                cigar = entry[6]
                bp = entry[7]
            out_line = "\t".join([chromosome, start, end, score, cigar, bp])
            output.write(out_line + "\n")


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
    bed_out,
    bed_in,
    ref,
    reads,
    genome,
    thread,
    out_dir,
    sample_prefix,
    min_mapq,
    min_af,
    method="max",
    slop=150,
):
    # # create masked genome (only allow bases around non-ref TEs)
    # expand_bed = out_dir + "/" + sample_prefix + ".nonref.expand.bed"
    # with open(expand_bed, "w") as output:
    #     subprocess.call(
    #         ["bedtools", "slop", "-i", bed_in, "-g", genome, "-b", str(slop)],
    #         stdout=output,
    #     )

    # complement_bed = out_dir + "/" + sample_prefix + ".nonref.complement.bed"
    # with open(complement_bed, "w") as output:
    #     subprocess.call(
    #         ["bedtools", "complement", "-i", expand_bed, "-g", genome], stdout=output
    #     )

    # masked_ref = ref + ".nonref.masked"
    # subprocess.call(
    #     ["bedtools", "maskfasta", "-fi", ref, "-bed", complement_bed, "-fo", masked_ref]
    # )
    # subprocess.call(["bwa", "index", masked_ref])
    # map raw reads to reference genome
    bam = out_dir + "/" + sample_prefix + ".nonref.bam"
    make_bam(fq=reads, ref=ref, thread=thread, bam=bam)

    # for each site, figure out AF
    samfile = pysam.AlignmentFile(bam, "rb")
    with open(bed_in, "r") as input, open(bed_out, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chromosome = entry[0]
            start = int(entry[1])
            end = int(entry[2])
            info = entry[3]
            score = entry[4]
            strand = entry[5]

            breakpoint_offset = 10

            tsd = int(info.split("|")[1])
            if tsd <= 0:
                SM_breakpoint = start
                MS_breakpoint = end
            else:
                SM_breakpoint = end
                MS_breakpoint = start

            n_cigar1 = count_nonref_read(
                samfile,
                chromosome,
                MS_breakpoint - breakpoint_offset,
                MS_breakpoint,
                "MS",
                min_mapq,
            )
            n_cigar2 = count_nonref_read(
                samfile,
                chromosome,
                SM_breakpoint,
                SM_breakpoint + breakpoint_offset,
                "SM",
                min_mapq,
            )

            ref_offset = 3
            n_ref = count_ref_read(
                samfile, chromosome, start - ref_offset, end + ref_offset, min_mapq
            )

            try:
                af1 = n_cigar1 / (n_cigar1 + n_ref)
            except ZeroDivisionError:
                af1 = None

            try:
                af2 = n_cigar2 / (n_cigar2 + n_ref)
            except ZeroDivisionError:
                af2 = None

            if method == "max":
                if af1 and af2:
                    af = round(max(af1, af2), 2)
                elif af1:
                    af = round(af1, 2)
                elif af2:
                    af = round(af2, 2)
                else:
                    af = "NA"
            else:
                if af1 and af2:
                    af = round(mean([af1, af2]), 2)
                elif af1:
                    af = round(af1, 2)
                elif af2:
                    af = round(af2, 2)
                else:
                    af = "NA"
            if af != "NA":
                if af >= min_af:
                    info_new = "|".join(
                        [info, str(af), str(n_cigar1), str(n_cigar2), str(n_ref)]
                    )
                    out_line = "\t".join(
                        [chromosome, str(start), str(end), info_new, score, strand]
                    )
                    output.write(out_line + "\n")


def count_nonref_read(samfile, chromosome, start, end, cigar_filter, min_mapq):
    n_read = 0
    breakpoint_error = 3
    for read in samfile.fetch(chromosome, start, end):
        if read.mapping_quality >= min_mapq:
            pattern, offset = parse_cigar(read.cigar)
            if pattern == cigar_filter:
                n_read = n_read + 1
            elif pattern == "SMS" and cigar_filter == "SM":
                read_breakpoint = read.reference_start
                if (
                    read_breakpoint >= start - breakpoint_error
                    and read_breakpoint <= start + breakpoint_error
                ):
                    n_read = n_read + 1
            elif pattern == "SMS" and cigar_filter == "MS":
                read_breakpoint = read.reference_start + offset
                if (
                    read_breakpoint >= end - breakpoint_error
                    and read_breakpoint <= end + breakpoint_error
                ):
                    n_read = n_read + 1
    return n_read


def count_ref_read(samfile, chromosome, start, end, min_mapq):
    n_read = 0
    for read in samfile.fetch(chromosome, start, end):
        if read.mapping_quality >= min_mapq:
            pattern, offset = parse_cigar(read.cigar)
            if pattern == "M":
                read_start = read.reference_start
                read_end = read_start + offset
                if read_start <= start and read_end >= end:
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


def refine_breakpoint(new_bed, old_bed, bam, cigar_type, min_mapq):
    samfile = pysam.AlignmentFile(bam, "rb")
    with open(old_bed, "r") as input, open(new_bed, "w") as output:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            chromosome = entry[0]
            cov_start = int(entry[1])
            cov_end = int(entry[2])
            breakpoints_sm = dict()
            breakpoints_ms = dict()
            cigars = {"SM": 0, "MS": 0}
            for read in samfile.fetch(chromosome, cov_start, cov_end):
                if read.mapping_quality >= min_mapq:
                    ref_start = read.reference_start
                    pattern, offset = parse_cigar(read.cigar)
                    if pattern == "SM":
                        cigars["SM"] = cigars["SM"] + 1
                        bp = ref_start
                        if bp in breakpoints_sm:
                            breakpoints_sm[bp] = breakpoints_sm[bp] + 1
                        else:
                            breakpoints_sm[bp] = 1
                    elif pattern == "MS":
                        cigars["MS"] = cigars["MS"] + 1
                        bp = ref_start + offset
                        if bp in breakpoints_ms:
                            breakpoints_ms[bp] = breakpoints_ms[bp] + 1
                        else:
                            breakpoints_ms[bp] = 1
                    else:
                        bp = None

            cigar_joint = get_keys(cigars)
            # breakpoint_joint = get_keys(breakpoints)

            if cigar_joint:
                if len(cigar_joint) > 1:
                    cigar_joint = "mix"
                else:
                    cigar_joint = cigar_joint[0]

            # if breakpoint_joint:
            #     breakpoint_joint = round(mean(breakpoint_joint))

            if cigar_joint == "MS":
                breakpoint_joint = round(mean(get_keys(breakpoints_ms)))
                refined_start = cov_start
                refined_end = breakpoint_joint
            elif cigar_joint == "SM":
                breakpoint_joint = round(mean(get_keys(breakpoints_sm)))
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
                        str(cigar_joint),
                        str(breakpoint_joint),
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
                    (int(entry[1]) - int(entry[7])) > 0
                    and (int(entry[2]) - int(entry[8])) > 0
                ) or (
                    (int(entry[7]) - int(entry[1])) > 0
                    and (int(entry[8]) - int(entry[2])) > 0
                ):  # get rid of entries if one is within another
                    if abs(int(entry[1]) - int(entry[8])) < abs(
                        int(entry[2]) - int(entry[7])
                    ):
                        if int(entry[1]) < int(entry[8]):
                            start = entry[1]
                            end = entry[8]
                        else:
                            start = entry[8]
                            end = entry[1]
                        strand = "-"
                        dist = int(entry[1]) - int(
                            entry[8]
                        )  # positive: gap; negative: overlap
                        if int(entry[1]) != int(entry[5]) or int(entry[8]) != int(
                            entry[11]
                        ):
                            ins_pass = False
                    else:
                        if int(entry[2]) < int(entry[7]):
                            start = entry[2]
                            end = entry[7]
                        else:
                            start = entry[7]
                            end = entry[2]
                        strand = "+"
                        dist = int(entry[7]) - int(
                            entry[2]
                        )  # positive: gap; negative: overlap
                        if int(entry[2]) != int(entry[5]) or int(entry[7]) != int(
                            entry[11]
                        ):
                            ins_pass = False
                    if dist < 0:
                        if -dist > tsd_max:
                            ins_pass = False
                    else:
                        if dist > gap_max:
                            ins_pass = False

                    if ins_pass:
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


def merge_bed(bed_in, bed_out, genome, filter_method="overlap"):
    # merge bed files from all families, check overlap, merge or remove entries if necessary
    bed_merge = bed_out + ".merge.tmp"
    with open(bed_merge, "w") as output:
        for bed in bed_in:
            if os.path.isfile(bed) and os.stat(bed).st_size != 0:
                with open(bed, "r") as input:
                    for line in input:
                        output.write(line)

    if get_lines(bed_merge) != 0:
        # sort bed files
        bed_sort = bed_out + ".merge.sort.tmp"
        with open(bed_sort, "w") as output:
            subprocess.call(
                ["bedtools", "sort", "-i", bed_merge, "-g", genome], stdout=output
            )
        os.remove(bed_merge)

        # remove entries if overlap with multiple families
        if filter_method == "overlap":
            bed_rm_overlap(bed_sort, bed_out)
        else:
            bed_rm_duplicate(bed_sort, bed_out)
        os.remove(bed_sort)

    else:
        os.rename(bed_merge, bed_out)


def parse_rm_out(rm_gff, gff3):
    with open(gff3, "w") as output, open(rm_gff, "r") as input:
        for line in input:
            if "RepeatMasker" in line:
                entry = line.replace("\n", "").split("\t")
                family = entry[8].split(" ")[1]
                family = re.sub('"Motif:', "", family)
                family = re.sub('"', "", family)
                out_line = "\t".join(
                    [
                        entry[0],
                        "RepeatMasker",
                        "dispersed_repeat",
                        entry[3],
                        entry[4],
                        entry[5],
                        entry[6],
                        entry[7],
                        "Target=" + family,
                    ]
                )
                output.write(out_line + "\n")


def gff3tobed(gff, bed):
    # check GFF3 format
    with open(gff, "r") as input:
        for line in input:
            if "#" not in line:
                if "Target=" not in line:
                    print(
                        "Incorrect GFF3 format, please check README for expected format, exiting..."
                    )
                    logging.exception(
                        "Incorrect GFF3 format, please check README for expected format, exiting..."
                    )
                    sys.exit(1)
                break
    with open(bed, "w") as output, open(gff, "r") as input:
        for line in input:
            if "#" not in line:
                entry = line.replace("\n", "").split("\t")
                info = entry[8].split(";")
                for item in info:
                    if "Target=" in item:
                        family = item.replace("Target=", "")
                out_line = "\t".join(
                    [entry[0], str(int(entry[3]) - 1), entry[4], family, ".", entry[6]]
                )
                output.write(out_line + "\n")


def get_ref(bed1, bed2, rm_bed, out_dir, family, window=50):
    # calculate clusters that jointly support ref TEs (all, norm) with a percentage
    ref_rm = out_dir + "/" + family + ".ref_rm.bed"
    family_ref_count = 0
    with open(ref_rm, "w") as output, open(rm_bed, "r") as input:
        for line in input:
            entry = line.replace("\n", "").split("\t")
            if entry[3] == family:
                family_ref_count = family_ref_count + 1
                output.write(line)

    if family_ref_count > 0:
        # calculate clusters that jointly support ref TEs (all, norm) with a percentage
        ref_sm = bed1.replace(".bed", ".ref.bed")
        if os.path.isfile(bed1):
            with open(ref_sm, "w") as output:
                subprocess.call(
                    [
                        "bedtools",
                        "window",
                        "-w",
                        str(window),
                        "-a",
                        ref_rm,
                        "-b",
                        bed1,
                        "-u",
                    ],
                    stdout=output,
                )
        ref_ms = bed2.replace(".bed", ".ref.bed")
        if os.path.isfile(bed2):
            with open(ref_ms, "w") as output:
                subprocess.call(
                    [
                        "bedtools",
                        "window",
                        "-w",
                        str(window),
                        "-a",
                        ref_rm,
                        "-b",
                        bed2,
                        "-u",
                    ],
                    stdout=output,
                )
        if os.path.isfile(ref_sm) and os.path.isfile(ref_ms):
            ref_both = out_dir + "/" + family + ".reference.bed"
            with open(ref_both, "w") as output:
                subprocess.call(
                    ["bedtools", "intersect", "-a", ref_sm, "-b", ref_ms, "-u"],
                    stdout=output,
                )
        if os.path.isfile(ref_sm):
            os.remove(ref_sm)
        if os.path.isfile(ref_ms):
            os.remove(ref_ms)
