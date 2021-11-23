#!/usr/bin/env python

import argparse
import os

# from genome_annotation import GenMakerCTL
from genome_annotation import GetLineage
from genome_annotation import exonerate_parser


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
    # jobfile.close()
    # os.system()
    # if busco_db_path == "":
    #     # print("wget %s%s" %(busco_web_path,buscodb.lower()))
    #     os.system("wget %s%s" % (busco_web_path, buscodb.lower()))
    # else:
    #     cpcmd = "cp %s%s ./" % (busco_db_path, buscodb.lower())
    #     os.system(cpcmd)


def run_snap():
    snap_exe = os.popen("which snap").read().replace("\n", "")
    snap_path = "/".join(snap_exe.split("/")[:-1]) + "/"
    perl_path = os.popen("which perl").read().replace("\n", "")
    # default snap_path set for genek user!
    os.system("mkdir snap")
    os.chdir("snap")
    # in path/to/genome/snap_rnd1
    os.system("mkdir training")
    os.system("cp ../exonerate/genome.fasta training/genome.dna")

    os.system("cp genome.fasta training")
    os.system("mkdir prediction")
    os.system("cp genome.fasta prediction")
    os.chdir("training")
    # os.system("maker2zff -n -d ../../maker_rnd1/genome.maker.output/genome_master_datastore_index.log")
    os.system(snap_path + "fathom genome.ann genome.dna -gene-stats ")
    os.system(snap_path + "fathom genome.ann genome.dna -validate")
    os.system(snap_path + "fathom genome.ann genome.dna -categorize 1000")
    os.system(snap_path + "fathom uni.ann uni.dna -export 1000 -plus")
    os.system("mkdir params")
    os.chdir("params")
    os.system(snap_path + "forge ../export.ann ../export.dna")
    os.chdir("..")
    os.system(
        perl_path + " " + snap_path + "hmm-assembler.pl genome.fasta params > genome.hmm")
    os.system("cp genome.hmm ../prediction")
    os.chdir("../prediction")
    os.system(
        snap_exe + " -aa genome.snap.protein.fasta -tx genome.snap.transcripts.fasta -name snap genome.hmm genome.fasta -gff > genome.snap.gff")
    return "This round of SNAP finished!"


def GetCDSgff():
    os.system("grep -v scoring genome.snap.gff > genome.snap.gff3")
    dummycheck()
    if os.path.exists('./dummy.txt'):
        os.system("rm dummy.txt")
        print("genome.snap.gffs saved!")
    os.system("sed 's/Einit/CDS/g; s/Esngl/CDS/g; s/Eterm/CDS/g ; s/Exon/CDS/g' genome.snap.gff > cds.gff")
    dummycheck()
    if os.path.exists('./dummy.txt'):
        os.system("rm dummy.txt")
        print("cds.gff saved!")
    path = os.getcwd()
    print(path + "/cds.gff")
    return path + "/cds.gff"


def Prep_antismash(genome, cdsgff):
    os.mkdir("antismash")
    os.chdir("antismash")
    # os.system("source activate antismash")
    os.system("cp %s ./genome.fasta" % (genome))
    os.system("cp %s ./" % (cdsgff))
    with open("smash.sh", 'w+') as smash:
        smash.write(
            "source activate antismash\nnohup antismash --taxon fungi --genefinding-tool none --genefinding-gff3 cds.gff genome.fasta &")
    smash.close()
    os.system("chmod +x smash.sh")
    cwd = os.getcwd()
    print("antismash script:" + cwd + "smash.sh")


def check_exe():
    pass


def annotate(spename, genomefilename, thread_num):
    # in path/to/genome/
    taxonclass = GetLineage.taxonkitWrapper(spename)
    os.system("mkdir exonerate")
    os.system("cp %s exonerate/genome.fasta" % (genomefilename))
    os.chdir("exonerate")
    # in path/to/genome/maker_rnd1
    exonerate_parser.run_exonerate(taxonclass, genomefilename, thread_num)

    os.chdir("../")
    # in path/to/genome/
    run_snap()


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    spename = args.name.replace("_", " ")
    taxonclass = GetLineage.taxonkitWrapper(spename)
    genomefilename = args.genome_file_name
    thread_num = str(args.thread_num)

    M_path = os.getcwd()

    annotate(spename, genomefilename, thread_num)
    gffpath = GetCDSgff()
    os.chdir(M_path)

    # Prep_antismash(M_path+"/"+genomefilename,gffpath)
