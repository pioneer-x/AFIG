#!/usr/bin/env python


def dummycheck():
    dummy = open("dummy.txt", 'w')
    dummy.write("It's dummy!")
    dummy.close()

def run_exonerate(taxonclass, genome, threads):
    import os
    import sys
    from joblib import Parallel, delayed
    script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    buscodb = os.path.join(script_path, "fungi_busco_sub" ,taxonclass + "_odb10.fasta")
    # jobfile = open("parallel_exon.sh", "w+")
    cmd_list = []
    for i in range(1, int(threads) + 1):
        cmd = "exonerate %s %s --querychunkid %s --querychunktotal %s --model protein2genome --showtargetgff yes --showvulgar no --showalignment no > %s.exon" %(buscodb, genome, str(i), str(threads), str(i))
        cmd_list.append(cmd)
        # jobfile.write((cmd))
    print("==" * 20)
    print("[INFO:] Running exonerate in protein2genome model.")
    Parallel(n_jobs=int(threads))(delayed(os.system)(cmd) for cmd in cmd_list)
    dummycheck()
    print("[INFO:] Done.")
    print("==" * 20)
    os.system("rm dummy.txt")
    os.system("cat *.exon > exon.out")
    exonerate_out = open("exon.out").read().split("#")
    return exonerate_out


def parse_exonerate(exonerate_out):
    zff_dict = {}
    gene_num = 0
    for i in range(15, len(exonerate_out)+1, 17):
        gene_num += 1
        gene_block = exonerate_out[i].split("\n")[1:-1]
        # print(exonerate_out[15])
        if len(gene_block) == 4:
            # Esngl
            for line in gene_block:
                line_list = line.split("\t")
                if line_list[2] == "exon":
                    start, end = line_list[3:5]
                if line_list[2] == "gene":
                    sequence = ">" + line_list[0]
                    model = line_list[8].split(";")[1].split(" ")[2] + "_" + str(i)
                    identity = float(line_list[8].split(";")[3].split(" ")[2])
            if sequence in zff_dict:
                zff_dict[sequence].append(["Esngl ", start, end, model, identity])
            else:
                zff_dict[sequence] = [["Esngl ", start, end, model, identity]]

            # zff_list.append(["Esngl", start, end, model])
        else:
            for line in gene_block:
                line_list = line.split("\t")
                if line_list[2] == "gene":
                    model = line_list[8].split(";")[1].split(" ")[2] + "_" + str(i)
                    identity = float(line_list[8].split(";")[3].split(" ")[2])
                    sequence = ">" + line_list[0]
                    gene_start, gene_end = line_list[3:5]
                if int(gene_start) < int(gene_end):
                    if line_list[2] == "exon":
                        start, end = line_list[3:5]
                        if start == gene_start:
                            # zff_list.append(["Einit", start, end, model])
                            if sequence in zff_dict:
                                zff_dict[sequence].append(["Einit ", start, end, model])
                            else:
                                zff_dict[sequence] = [["Einit ", start, end, model]]
                        if end == gene_end:
                            # zff_list.append(["Eterm", start, end, model])
                            if sequence in zff_dict:
                                zff_dict[sequence].append(["Eterm ", start, end, model])
                            else:
                                zff_dict[sequence] = [["Eterm ", start, end, model]]
                        if start != gene_start and end != gene_end:
                            # zff_list.append(["Exon", start, end, model])
                            if sequence in zff_dict:
                                zff_dict[sequence].append(["Exon ", start, end, model])
                            else:
                                zff_dict[sequence] = [["Exon ", start, end, model]]
                if int(gene_start) > int(gene_end):
                    if line_list[2] == "exon":
                        start, end = line_list[3:5]
                        if start == gene_start:
                            # zff_list.append(["Einit", start, end, model])
                            if sequence in zff_dict:
                                zff_dict[sequence].append(["Eterm ", start, end, model])
                            else:
                                zff_dict[sequence] = [["Eterm ", start, end, model]]
                        if end == gene_end:
                            # zff_list.append(["Eterm", start, end, model])
                            if sequence in zff_dict:
                                zff_dict[sequence].append(["Einit ", start, end, model])
                            else:
                                zff_dict[sequence] = [["Einit ", start, end, model]]
                        if start != gene_start and end != gene_end:
                            # zff_list.append(["Exon", start, end, model])
                            if sequence in zff_dict:
                                zff_dict[sequence].append(["Exon ", start, end, model])
                            else:
                                zff_dict[sequence] = [["Exon ", start, end, model]]
    return zff_dict

def write_zff(zff_dict):
    with open("genome.ann", 'w+') as zff_out:
        for sequence in zff_dict:
            zff_out.write(sequence + "\n")
            for line in zff_dict[sequence]:
                zff_out.write(" ".join(line) + "\n")
        zff_out.close()

if __name__ == '__main__':
    import sys
    from joblib import Parallel, delayed
    # print(sys.argv[1])
    # run_exonerate(taxonclass, genome, threads)
    zff_dict = parse_exonerate(open(sys.argv[1]).read().split("#"))
    write_zff(zff_dict)
