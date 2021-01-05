#!/usr/bin/env python3

import argparse
import sys
import os
import time
import logging
import subprocess
import shutil
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
This script predicts non-reference TE insertions using single-end short read data
"""


def get_args():
    parser = argparse.ArgumentParser(
        description="Script to detect non-reference TEs from single end short read data"
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("required arguments")

    ## required ##
    required.add_argument(
        "-f",
        "--read",
        type=str,
        help="raw reads in fastq or fastq.gz format",
        required=True,
    )
    required.add_argument(
        "-l", "--library", type=str, help="TE concensus sequence", required=True
    )
    required.add_argument(
        "-r", "--reference", type=str, help="reference genome", required=True
    )
    ## optional
    optional.add_argument(
        "-n", "--region", type=str, help="region to filter", required=False
    )
    optional.add_argument(
        "-w",
        "--window",
        type=int,
        help="merge window for identifying TE clusters (default = 100bp) ",
        required=False,
    )
    optional.add_argument(
        "--experiment",
        action="store_true",
        help="If provided then reads will be mapped to masked augmented reference in the first step (by default reads will be mapped to TE library)",
        required=False,
    )
    optional.add_argument(
        "--tsd_max", type=int, help="maximum TSD size (default = 20) ", required=False
    )
    optional.add_argument(
        "--gap_max", type=int, help="maximum gap size (default = 0) ", required=False
    )
    optional.add_argument(
        "-m", "--mapper", type=str, help="mapper (default = bwa) ", required=False
    )
    optional.add_argument(
        "-t", "--thread", type=int, help="thread (default = 1) ", required=False
    )
    optional.add_argument(
        "-o", "--out", type=str, help="output dir (default = '.') ", required=False
    )
    optional.add_argument(
        "-k",
        "--keep_files",
        action="store_true",
        help="If provided then all intermediate files will be kept (default: remove intermediate files)",
        required=False,
    )
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # checks if in files exist
    try:
        test = open(args.read, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.read)
        sys.exit(1)

    try:
        test = open(args.reference, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.reference)
        sys.exit(1)

    try:
        test = open(args.library, "r")
    except Exception as e:
        print(e)
        logging.exception("Can not open input file: " + args.library)
        sys.exit(1)

    # sets up out dir variable
    if args.out is None:
        args.out = "."
    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        os.mkdir(args.out)

    if args.mapper is None:
        args.mapper = "bwa"

    if args.tsd_max is None:
        args.tsd_max = 20

    if args.gap_max is None:
        args.gap_max = 0

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
    start_time_all = time.time()

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

    # prepare modified reference genome
    if args.experiment:
        augment = True
    else:
        augment = False
    rm_dir = os.path.join(tmp_dir, "repeatmask")
    mkdir(rm_dir)
    ref_modified, te_gff = repeatmask(
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
    start_time_align = time.time()
    if not args.experiment:
        # align reads to TE library (single end mode)
        bam = tmp_dir + "/" + sample_name + ".bam"
        make_bam(fastq, library, str(args.thread), bam, args.mapper)
        os.remove(fastq)
    else:
        # align reads to masked augmented reference (single end mode)
        bam = tmp_dir + "/" + sample_name + ".bam"
        make_bam(fastq, ref_modified, str(args.thread), bam, args.mapper)
        os.remove(fastq)
    proc_time_align = time.time() - start_time_align
    logging.info("Alignment finished in " + format_time(proc_time_align))

    # use samtools to separate bam into family bam (single thread)
    print("Detection by family...")
    logging.info("Start insertion candidate search...")
    start_time_candidate = time.time()
    family_dir = os.path.join(tmp_dir, "family-specific")
    mkdir(family_dir)
    family_arguments = []
    for family in families:
        argument = [
            family,
            bam,
            ref_modified,
            te_gff,
            family_dir,
            args.mapper,
            contigs,
            args.experiment,
            args.tsd_max,
            args.gap_max,
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
    proc_time_candidate = time.time() - start_time_candidate
    logging.info(
        "Insertion candidate search finished in " + format_time(proc_time_candidate)
    )

    # merge non-ref bed files
    final_bed = args.out + "/" + sample_name + ".nonref.bed"
    pattern = "/*/*.nonref.bed"
    bed_files = glob(family_dir + pattern, recursive=True)
    merge_bed(bed_in=bed_files, bed_out=final_bed)

    # merge ref bed files
    if te_gff is not None:
        final_bed = args.out + "/" + sample_name + ".ref.bed"
        pattern = "/*/*.ref.bed"
        bed_files = glob(family_dir + pattern, recursive=True)
        merge_bed(bed_in=bed_files, bed_out=final_bed)

    # clean tmp files
    if not args.keep_files:
        shutil.rmtree(tmp_dir)

    proc_time_all = time.time() - start_time_all
    print("ngs_te_mapper finished!")
    logging.info("ngs_te_mapper finished in " + format_time(proc_time_all))


main()