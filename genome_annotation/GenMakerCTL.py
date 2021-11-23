#!/usr/bin/python


# To generate maker opt file for protein2genome annotation

import argparse
import os
import sys

busco_db_path = os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])), "fungi_busco_sub/")
# Server was down!
busco_web_path = "http://202.119.249.49/busco/"


def get_parser():
    parser = argparse.ArgumentParser(description='Get a busco taxon class and generate maker_opt file')
    # parser.parse_args()
    parser.add_argument("-g", "--genome_file_name", default="genome.fasta")
    parser.add_argument("-t", "--taxon_name", default="fungi")
    parser.add_argument("-p", "--thread_num", default=4)

    return parser


def GenMaker_opt(taxonclass, threads):
    buscodb = taxonclass + "_odb10.fasta"
    if busco_db_path == "":
        # print("wget %s%s" %(busco_web_path,buscodb.lower()))
        os.system("wget %s%s" % (busco_web_path, buscodb.lower()))
    else:
        cpcmd = "cp %s%s ./" % (busco_db_path, buscodb.lower())
        print(cpcmd)
        os.system(cpcmd)

    os.system("maker -CTL")
    sedcmd = "sed -i 's/^genome=/genome=genome.fasta/g;s/^repeat_protein=\/pub\/software\/maker-3.01.03\/data\/te_proteins.fasta/repeat_protein= /g;s/^model_org=all/model_org=0/g;s/^softmask=1/softmask=0/g;s/^cpus=1/cpus=%s/g;s/^protein2genome=0/protein2genome=1/g;s/^protein= /protein=%s/g' maker_opts.ctl" % (
    threads, buscodb.lower())
    os.system(sedcmd)
    return "Modiffied maker_opt file successfully!"


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    taxonclass = args.taxon_name
    genomefilename = args.genome_file_name
    os.system("cp %s genome.fasta" % (genomefilename))
    thread_num = str(args.thread_num)
    # print(taxonclass,genomefilename)
    GenMaker_opt(taxonclass, thread_num)


