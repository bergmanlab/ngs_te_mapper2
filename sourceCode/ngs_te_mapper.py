#!/usr/bin/env python3

import argparse
import sys
import os
import time
import logging
import subprocess
from glob import glob
from multiprocessing import Pool
from utility import (
    parse_input,
    repeatmask,
    make_bam,
    get_family_bed,
    merge_bed,
    mkdir,
    format_time,
)

"""
Author: Shunhua Han <shhan@uga.edu>
this script predicts non-reference TE insertions using single-end short read data
"""


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to detect non-reference TEs from single end short read data"
    )
    ## required ##
    parser.add_argument(
        "-f",
        "--read",
        type=str,
        help="raw reads in fastq or fastq.gz format",
        required=True,
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

    # logging config
    formatstr = "%(asctime)s: %(levelname)s: %(message)s"
    datestr = "%m/%d/%Y %H:%M:%S"
    logging.basicConfig(
        level=logging.DEBUG,
        filename=os.path.join(args.out, "ngs_te_mapper.log"),
        filemode="w",
        format=formatstr,
        datefmt=datestr,
    )
    logging.info("CMD: " + " ".join(sys.argv))
    start_time = time.time()

    # create directory for intermediate files
    tmp_dir = os.path.join(args.out, "intermediate_files")
    mkdir(tmp_dir)

    # Parse input
    sample_name = os.path.basename(args.read).replace(".fastq.gz", "")
    fastq, library, ref = parse_input(
        input_read=args.read,
        input_library=args.library,
        input_reference=args.reference,
        out_dir=tmp_dir,
    )
    print(fastq)
    print(library)
    print(ref)

    # prepare modified reference genome
    if args.ngs_te_mapper:
        augment = False
    else:
        augment = True
    rm_dir = os.path.join(tmp_dir, "repeatmask")
    mkdir(rm_dir)
    ref_modified = rm_dir + "/" + "dm6.fasta.masked"
    ref_modified = repeatmask(
        ref=ref,
        library=library,
        outdir=rm_dir,
        thread=args.thread,
        augment=False,
    )
    if args.mapper == "bwa":
        subprocess.call(["bwa", "index", ref_modified])
        subprocess.call(["bwa", "index", library])

    # get all TE families
    families = []
    with open(library, "r") as input:
        for line in input:
            if ">" in line:
                family = line.replace("\n", "").replace(">", "")
                families.append(family)

    # get all contigs from reference
    contigs = []
    with open(ref, "r") as input:
        for line in input:
            if ">" in line:
                contig = line.replace("\n", "").replace(">", "")
                contigs.append(contig)

    contigs = " ".join(contigs)

    # step one: align read to TE library or masked augmented ref
    print("Align reads to TE library..")
    logging.info("Start alignment...")
    start_time = time.time()
    if args.ngs_te_mapper:
        # align reads to TE library (single end mode)
        bam = tmp_dir + "/" + sample_name + ".bam"
        make_bam(fastq, library, str(args.thread), bam, args.mapper)
        os.remove(fastq)
    else:
        # align reads to masked augmented reference (single end mode)
        bam = tmp_dir + "/" + sample_name + ".bam"
        make_bam(fastq, ref_modified, str(args.thread), bam, args.mapper)
        os.remove(fastq)
    proc_time = time.time() - start_time
    logging.info("Alignment finished in " + format_time(proc_time))

    # use samtools to separate bam into family bam (single thread)
    print("Detection by family...")
    family_dir = os.path.join(tmp_dir, "family-specific")
    mkdir(family_dir)
    family_arguments = []
    for family in families:
        argument = [
            family,
            bam,
            ref_modified,
            family_dir,
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
    bed_files = glob(family_dir + pattern, recursive=True)
    merge_bed(bed_in=bed_files, bed_out=final_bed)

    proc_time = time.time() - start_time
    print("ngs_te_mapper finished!")
    logging.info("ngs_te_mapper finished in " + format_time(proc_time))


main()