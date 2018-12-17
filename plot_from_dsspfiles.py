#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import glob
import math
import matplotlib.pyplot as plt

from collections import Counter, defaultdict


CysLinkChange = "abcdefghijklmnopqrstuvwxyz"
ResCode = {
    'ALA':"A",
    'ASN':"N",
    'ASP':"D",
    'ARG':"R",
    'CYS':"C",
    'GLN':"Q",
    'GLU':"E",
    'GLY':"G",
    'HIS':"H",
    'ILE':"I",
    'LEU':"L",
    'LYS':"K",
    'MET':"M",
    'PRO':"P",
    'PHE':"F",
    'SER':"S",
    'THR':"T",
    'TRP':"W",
    'TYR':"Y",
    'VAL':"V",
    "A":'ALA',
    "N":'ASN',
    "D":'ASP',
    "R":'ARG',
    "C":'CYS',
    "Q":'GLN',
    "E":'GLU',
    "G":'GLY',
    "H":'HIS',
    "I":'ILE',
    "L":'LEU',
    "K":'LYS',
    "M":'MET',
    "P":'PRO',
    "F":'PHE',
    "S":'SER',
    "T":'THR',
    "W":'TRP',
    "Y":'TYR',
    "V":'VAL'
}


def read_chains(n_dssp):
    print(n_dssp)
    chains = defaultdict(list)
    with open(n_dssp, "r") as dsspfile:
        RL = dsspfile.readline()
        # skip dssp file header
        while not "#  RESIDUE AA STRUCTURE" in RL:
            RL = dsspfile.readline()
        RL = dsspfile.readline()

        while RL:
            try:
                # #         numdssp   numpdb   Chain     AA       S        phi             psy
                # #        [RL[0:5], RL[5:10], RL[11], RL[13], RL[16], RL[103:109], RL[109:115]]
                readAmin = [int(RL[0:5]), int(RL[5:10]), RL[11], RL[13], RL[16], RL[103:109], RL[109:115]]
            except ValueError:
                print("!!!!Err Decoding dssp record", RL)
            else:
                chains[readAmin[2]].append(readAmin)
            RL = dsspfile.readline()
    return chains


def append_tetra_oligs(tetra_oligs, chains):
    for l_chain in chains.keys():
        for iAmin in range(len(chains[l_chain])-3):
            subchain = chains[l_chain][iAmin:iAmin+4]
            # check if elements are aligned in order (ex. 1-2-3-4, not 1-5-8-9)
            if [x[1]-chains[l_chain][iAmin][1] for x in subchain] == [0,1,2,3]:
                aminseq = ""
                for AminCheck in [x[3] for x in subchain]:
                    if AminCheck in ResCode.keys():
                        aminseq += AminCheck
                    if AminCheck in CysLinkChange:
                        aminseq += "C"
                if len(aminseq) == 4:
                    tetra_oligs[aminseq] += 1


def write_tetras(tetra_oligs):
    with open("tet_oligs", "w") as tetra_oligs_file:
        sortoligs = sorted(tetra_oligs.keys())
        for oliga in sortoligs:
            tetra_oligs_file.write(oliga + ":" + str(tetra_oligs[oliga]) + "\n")


if __name__ == "__main__":
    tetra_oligs = defaultdict(int)

    for n_dssp in glob.glob("./*.dssp"): #".../*.dssp"):  #path to the directory in second quotes
        chains = read_chains(n_dssp)
        append_tetra_oligs(tetra_oligs, chains)

    write_tetras(tetra_oligs)

    Cnt = Counter(list(tetra_oligs.values()))
    X = sorted(Cnt)
    Y = [Cnt[x] for x in X]
    logY = [math.log(y) for y in Y]
    plt.plot(X, logY)
    plt.show()
