#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    if(isfile(amplicon_file)):
        with gzip.open(amplicon_file, 'r') as f:
            file_content = f.read().decode('utf-8')
            lines = iter(str(file_content).split('\n'))
            liste=[]
            seq=''
            for line in lines:
                if len(line) != 0:
                    if line[0]=='>':
                        if seq != '':
                            liste.append(seq)
                        seq=''
                    else:
                        seq=seq+line
            liste.append(seq)
            for seq in iter(liste):
                if len(seq)>=minseqlen:
                    yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dic = {}
    for gen in read_fasta(amplicon_file,minseqlen):
        if gen in dic:
            dic[gen]+=1
        else:
            dic[gen]=1

    for keys,value in sorted(dic.items(),key=lambda x: x[1],reverse=True):
        if value >= mincount:
            yield [keys,value]


def get_chunks(sequence, chunk_size):
    L=[]
    j=0
    temp=""
    for i in range(len(sequence)):
        temp+=sequence[i]
        if (j==chunk_size-1):
            L.append(temp)
            temp=""
            j=0
        else:
            j+=1
    try:
        len(L)>=4
    except:
        raise ValueError
    else:
        return L


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence) - kmer_size +1):
        yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
    c=0
    for i in range(len(alignment_list[0])):
        if (alignment_list[0][i]==alignment_list[1][i]):
            c+=1
    return c/len(alignment_list[0])*100

def get_unique_kmer(kmer_dict ,sequence,id_seq,kmer_size):
    kmer_gen = cut_kmer(sequence,kmer_size)
    for kmer in kmer_gen:
        if (kmer_dict.get(kmer)==None):
            kmer_dict[kmer]=kmer_dict.get(kmer,[id_seq])
        elif (id_seq not in kmer_dict.get(kmer) ):
            kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    cnt = Counter()
    kmer_gen = cut_kmer( sequence, kmer_size)
    for kmer in kmer_gen:
        if kmer in kmer_dict:
            #print(kmer_dict[kmer])
            for i in kmer_dict[kmer]:
                cnt[i] += 1
    mostcommons = cnt.most_common(8)
    L=[]
    for i in mostcommons:
        L.append(i[0])
    return L

def std(data):
    return statistics.stdev(data)

def detect_chimera(perc_identity_matrix):
    L=[]
    dif_similarites = 0
    for l in perc_identity_matrix :
        L.append(std(l))
        if (l[0] != perc_identity_matrix[0][0] or l[1] != perc_identity_matrix[0][1]):
            dif_similarites +=1
    ecart_type_moyen=(statistics.mean(L))
    return (ecart_type_moyen >5 and dif_similarites >= 2)


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    gen_seq = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    kmer_dict ={}
    id_seq = 0
    sequence_list = []
    
    for sequence in gen_seq:
        perc_identity_matrix = []
        sequence_list.append(sequence[0])
        segments = get_chunks(sequence[0],chunk_size)

        
        
        #division de chaque sequence candidate en 4 segments de longueur chunk_size
        for seg in segments:

            kmer_dict = get_unique_kmer(kmer_dict,seg,id_seq,kmer_size)
        mates = search_mates(kmer_dict, sequence[0], kmer_size)

        if len(mates)<=2 :
            yield sequence
        else:
            parents=[get_chunks(sequence_list[mates[0]],chunk_size),get_chunks(sequence_list[mates[1]],chunk_size)]
            
                
            for j,chunk in enumerate(segments):
                L=[]
                for i in range(len(parents)):

                    text=nw.global_align(chunk, parents[i][j])
                    r = get_identity(text)

                    L.append(r)
                perc_identity_matrix.append(L)

            if not detect_chimera(perc_identity_matrix):

                yield sequence
        #faire la matrice pour detect chimera
        
        id_seq +=1

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    OTU = []
    lst = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    L=[]
    first=True
    for seq  in lst:
        if (first):
            OTU.append(seq)
            first = False
        else:
            for j in range(len(L)):
                if (get_identity(nw.global_align(seq[0],L[j][0])) <= 97) :
                    OTU.append(seq)
        L.append(seq)
    return OTU


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file, "w") as f:
        for k, _ in enumerate(OTU_list):
            f.write(">OTU_" + str(k + 1) + " occurrence:"+ str(OTU_list[k][1]) + "\n")
            f.write(fill(str(OTU_list[k][0])))
            f.write("\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()