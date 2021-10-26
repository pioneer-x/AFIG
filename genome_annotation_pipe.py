#!/usr/bin/env python

import argparse
import os

from genome_annotation import GenMakerCTL
from genome_annotation import GetLineage


def get_parser():
    parser = argparse.ArgumentParser(description='Automate, RNA-seq free, fungi genome annotation pipeline.',
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog=(
                                         "demo: python3 genome_annotation_pipe.py -n Apiotrichum_mycotoxinovorans -g Am_genome.fasta -p 16"))
    # parser.parse_args()
    parser.add_argument("-n", "--name")
    parser.add_argument("-g", "--genome_file_name", default="genome.fasta")
    parser.add_argument("-p", "--thread_num", default=16)
    return parser


def dummycheck():
    dummy = open("dummy.txt", 'w')
    dummy.write("It's dummy!")
    dummy.close()
    # if os.path.exists('./dummy.txt'):
    # return "OK"

def run_snap():
    os.system("mkdir snap_rnd1")
    os.chdir("snap_rnd1")
    #in path/to/genome/snap_rnd1
    os.system("cp ../maker_rnd1/genome.fasta ./")
    os.system("mkdir training")
    os.system("cp genome.fasta training")
    os.system("mkdir prediction")
    os.system("cp genome.fasta prediction")
    os.chdir("training")
    os.system("maker2zff -n -d ../../maker_rnd1/genome.maker.output/genome_master_datastore_index.log")
    os.system("/pub/software/maker-2.31.10/exe/snap/fathom genome.ann genome.dna -gene-stats ")
    os.system("/pub/software/maker-2.31.10/exe/snap/fathom genome.ann genome.dna -validate")
    os.system("/pub/software/maker-2.31.10/exe/snap/fathom genome.ann genome.dna -categorize 1000")
    os.system("/pub/software/maker-2.31.10/exe/snap/fathom uni.ann uni.dna -export 1000 -plus")
    os.system("mkdir params")
    os.chdir("params")
    os.system("/pub/software/maker-2.31.10/exe/snap/forge ../export.ann ../export.dna")
    os.chdir("..")
    os.system("/pub/perl/bin/perl pub/software/maker-2.31.10/exe/snap/hmm-assembler.pl genome.fasta params > genome.hmm")
    os.system("cp genome.hmm ../prediction")
    os.chdir("../prediction")
    os.system(
        "/pub/software/maker-2.31.10/exe/snap/snap -aa genome.snap.rnd1.protein.fasta -tx genome.snap.rnd1.transcripts.fasta -name snap genome.hmm genome.fasta -gff > genome.snap.rnd1.gff")
    return "This round of SNAP finished!"

def GetCDSgff():
    os.system("grep -v scoring genome.snap.rnd1.gff > genome.snap.rnd1.gff3")
    dummycheck()
    if os.path.exists('./dummy.txt'):
        os.system("rm dummy.txt")
        print("genome.snap.rnd1.gffs saved!")
    os.system("sed 's/Einit/CDS/g; s/Esngl/CDS/g; s/Eterm/CDS/g ; s/Exon/CDS/g' genome.snap.rnd1.gff > cds.gff")
    dummycheck()
    if os.path.exists('./dummy.txt'):
        os.system("rm dummy.txt")
        print("cds.gff saved!")
    path = os.getcwd()
    print(path+"/cds.gff")
    return path+"/cds.gff"

def Prep_antismash(genome,cdsgff):
    os.mkdir("antismash")
    os.chdir("antismash")
    # os.system("source activate antismash")
    os.system("cp %s ./genome.fasta" %(genome))
    os.system("cp %s ./" %(cdsgff))
    with open("smash.sh", 'w+') as smash:
        smash.write(
            "source activate antismash\nnohup antismash --taxon fungi --genefinding-tool none --genefinding-gff3 cds.gff genome.fasta &")
    smash.close()
    os.system("chmod +x smash.sh")
    cwd = os.getcwd()
    print("antismash script:" + cwd + "smash.sh")


def annotate(spename,genomefilename,thread_num):
    #in path/to/genome/
    taxonclass = GetLineage.taxonkitWrapper(spename)
    os.system("mkdir maker_rnd1")
    os.system("cp %s maker_rnd1/genome.fasta" % (genomefilename))
    os.chdir("maker_rnd1")
    #in path/to/genome/maker_rnd1
    GenMakerCTL.GenMaker_opt(taxonclass, thread_num)
    # MPI not supportted
    # os.system("mpiexec -n %s maker" % (thread_num))
    os.system("maker")
    os.chdir("genome.maker.output")
    #in path/to/genome/maker_rnd1/genome.maker.output/
    os.system("gff3_merge -s -d genome_master_datastore_index.log > genome_rnd1.busco.maker.gff")
    os.system("fasta_merge -d genome_master_datastore_index.log")
    os.chdir("../../")
    #in path/to/genome/
    run_snap()


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    spename = args.name.replace("_", " ")
    taxonclass = GetLineage.taxonkitWrapper(spename)
    genomefilename = args.genome_file_name
    thread_num = str(args.thread_num)

    M_path = os.getcwd()

    annotate(spename,genomefilename,thread_num)
    gffpath = GetCDSgff()
    os.chdir(M_path)

    Prep_antismash(M_path+"/"+genomefilename,gffpath)




