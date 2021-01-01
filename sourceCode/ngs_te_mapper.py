#!/usr/bin/env python3

import argparse
import sys
import os
import subprocess
from glob import glob
from multiprocessing import Process, Pool


# from Bio import SeqIO
from utility import (
    make_bam,
    get_family_bed,
    merge_bed,
    mkdir,
)

"""
this script takes single end illumina data and output reference/non-reference TEs for input TE families
"""


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to detect non-reference TEs from single end short read data"
    )
    ## required ##
    parser.add_argument(
        "-f", "--reads", type=str, help="fastq.gz file or file list", required=True
    )
    parser.add_argument(
        "-l", "--library", type=str, help="TE concensus sequence", required=True
    )
    parser.add_argument(
        "-r", "--reference", type=str, help="reference genome", required=True
    )
    ## optional
    parser.add_argument(
        "-n", "--region", type=str, help="region to filter", required=False
    )
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        help="merge window for identifying TE clusters (default = 100bp) ",
        required=False,
    )
    parser.add_argument(
        "--ngs_te_mapper",
        action="store_true",
        help="If provided then reads will be mapped to TE in the first step (like in ngs_te_mapper)",
        required=False,
    )
    parser.add_argument(
        "-m", "--mapper", type=str, help="mapper (default = bwa) ", required=False
    )
    parser.add_argument(
        "-t", "--thread", type=int, help="thread (default = 1) ", required=False
    )
    parser.add_argument(
        "-o", "--out", type=str, help="output dir (default = '.') ", required=False
    )
    args = parser.parse_args()

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.mapper is None:
        args.mapper = "bwa"

    if args.thread is None:
        args.thread = 1

    if args.window is None:
        args.window = 100

    return args


def main():
    args = get_args()

    # reads list
    raw_reads = args.reads.replace(" ", "").split(",")
    sample_name = os.path.basename(raw_reads[0]).replace(".fastq.gz", "")

    # unzip and merge input files, if muliple inputs were provided
    fastq = args.out + "/" + sample_name + ".fastq"
    with open(fastq, "w") as output:
        for read in raw_reads:
            subprocess.call(["gunzip", "-c", read], stdout=output)

    # ## create masked augmented reference genome
    # ref_new = "/scratch/sh60271/ngs2/data/dm6.cov.fasta.masked.aug"
    # ref_rm = "/scratch/sh60271/ngs2/data/dm6.rm.fasta"


    # prepare modified reference genome
    if args.ngs_te_mapper:
        augment = False
    else:
        augment = True
    rm_dir = os.path.join(args.out, "repeatmask")
    mkdir(rm_dir)
    ref_modified = repeatmask(
        ref=args.reference,
        library=args.library,
        outdir=rm_dir,
        thread=args.thread,
        augment=False,
    )
    if args.mapper == "bwa":
        subprocess.call(["bwa", "index", ref_modified])

    # process families
    families = []
    with open(args.library, "r") as input:
        for line in input:
            if ">" in line:
                family = line.replace("\n", "").replace(">", "")
                families.append(family)

    # process reference
    contigs = []
    with open(args.reference, "r") as input:
        for line in input:
            if ">" in line:
                contig = line.replace("\n", "").replace(">", "")
                contigs.append(contig)

    contigs = " ".join(contigs)

    # step one: read alignment to TE library or masked augmented ref
    print("align reads to TE library..")
    if args.ngs_te_mapper:
        # align reads to TE library (single end mode)
        bam = args.out + "/" + sample_name + ".bam"
        make_bam(fastq, args.library, str(args.thread), bam, args.mapper)
        os.remove(fastq)
    else:
        # align reads to masked augmented reference (single end mode)
        bam = args.out + "/" + sample_name + ".bam"
        make_bam(fastq, ref_new, str(args.thread), bam, args.mapper)
        os.remove(fastq)

    # use samtools to separate bam into family bam (single thread)
    print("detection by family...")
    family_arguments = []
    for family in families:
        argument = [
            family,
            bam,
            ref_rm,
            args.out,
            args.mapper,
            contigs,
            args.ngs_te_mapper,
        ]
        family_arguments.append(argument)
    try:
        pool = Pool(processes=args.thread)
        pool.map(get_family_bed, family_arguments)
        pool.close()
        pool.join()
    except Exception as e:
        print(e)
        print("Family detection failed, exiting...")
        sys.exit(1)

    # merge non ref bed files
    final_bed = args.out + "/" + sample_name + ".nonref.bed"
    pattern = "/*/*nonref.bed"
    bed_files = glob(args.out + pattern, recursive=True)
    merge_bed(bed_in=bed_files, bed_out=final_bed)


main()