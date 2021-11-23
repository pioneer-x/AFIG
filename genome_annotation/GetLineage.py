#!/usr/bin/python

# To get the lineage by name using taxonkit for busco chosing

import argparse

from ete3 import NCBITaxa


def get_parser():
    parser = argparse.ArgumentParser(description='To get the lineage by name using taxonkit and chose busco dataset')
    # parser.parse_args()
    parser.add_argument("-n", "--name")
    return parser


def taxonWrapper(taxon_name):

    busco_list = [
        "agaricales",
        "agaricomycetes",
        "ascomycota",
        "basidiomycota",
        "boletales",
        "capnodiales",
        "chaetothyriales",
        "dothideomycetes",
        "eurotiales",
        "eurotiomycetes",
        "fungi",
        "glomerellales",
        "helotiales",
        "hypocreales",
        "leotiomycetes",
        "microsporidia",
        "mucorales",
        "mucoromycota",
        "onygenales",
        "pleosporales",
        "polyporales",
        "saccharomycetes",
        "sordariomycetes",
        "tremellomycetes",
    ]
    ncbi = NCBITaxa()
    tax_id = ncbi.get_name_translator([taxon_name])[taxon_name][0]
    lineage = ncbi.get_lineage(tax_id)
    lineage_list = [taxname.lower() for taxname in ncbi.get_taxid_translator(lineage).values()]
    for taxonclass in lineage_list[::-1]:
        if taxonclass in busco_list:
            print("==" * 20)
            print("[INFO:] Found %s(taxid: %s) within %s lineage." %(taxon_name, str(tax_id), taxonclass))
            print("==" * 20)
            return taxonclass


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    spename = args.name.replace("_", " ")
    taxonname = taxonWrapper(spename)
    print(taxonname)
