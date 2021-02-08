#! /usr/bin/env python3
# Author: Kristopher Kieft, UW-Madison
# Collaborator: Zhichao Zhou, UW-Madison


# IMPORTANT
# All proteins must have a prefix respective to their location/sample origin

# Works best for inputs of at least 5 locations/samples

# Core proteins: clusters containing at least 1 protein from at least "-l" different locations. This will result in a minimum cluster size of "-l" proteins if each protein comes from a unique location
# Unique proteins: clusters containing at least "-c" proteins from a single location, must be unique phages up to "-c" minimum (can't have 5 proteins from 1 phage count)
# Flexible proteins: clusters containing at least "-c" proteins from at least 2 locations but no more than "-l" minus one locations, must be unique phages up to "-c" minimum (can't have 5 proteins from 1 phage count)


try:
    import sys
    import subprocess
    from Bio.SeqIO.FastaIO import SimpleFastaParser as SFP
    import argparse
    import time
    import datetime
    from datetime import date
    import logging
    import uuid
    import os
except Exception as e:
    sys.stderr.write("\nError: please verify dependancy imports are installed and up to date:\n\n")
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

start = time.time()
start_time = str(datetime.datetime.now().time()).rsplit(".",1)[0]
u = str(uuid.uuid1()).split("-")[0]

# PhageCore
PhageCore = argparse.ArgumentParser(description='identifies core, flexible and unique protein clusters from mmseqs2 results')
PhageCore.add_argument('--version', action='version', version='PhageCore v0.0.0')
PhageCore.add_argument('-i', type=str, required=True, nargs=1, help='input mmseqs2 "all_seqs.fasta" result [required]')
PhageCore.add_argument('-o', type=str, nargs=1, required=True, help="name of output folder (will create if doesn't exist) [required]")
PhageCore.add_argument('-l', type=str, nargs=1, default = ['9'], help='minimum number of samples/locations for "core" proteins (integer) [default = 9]')
PhageCore.add_argument('-c', type=str, nargs=1, default = ['5'], help='minimum cluster size for "flexible" and "unique" proteins (integer) [default = 5]')
PhageCore.add_argument('-p', type=str, nargs=1, default = ['~~'], help='string deliminating location/sample prefix. If prefix contains a dash (-) then use quotes [default = ~~]')
PhageCore.add_argument('-g', type=str, nargs=1, default = ['2'], help='minimum number of unique genomes per cluster [default = 2]')
PhageCore.add_argument('-m', type=str, nargs=1, default = ['1'], help='minimum number of proteins per cluster per location for flexible and core [default = 1]')
PhageCore.add_argument('--force', action='store_true', help='force write to existing output folder (may mess up clusters) [default = off]')
args = PhageCore.parse_args()
#
mmseqs = str(args.i[0])
folder = str(args.o[0])
if folder[-1] != "/":
    folder += "/"
if os.path.exists(str(folder)) and args.force == False:
    print("\nError: please remove existing output folder or set new output folder name. Exiting.\n")
    exit(1)
if os.path.exists(str(folder)) and args.force == True:
    print("\nWriting to existing output folder. Caution: will mess up results if clusters already exist.\n")
if not os.path.exists(str(folder)):
    subprocess.run('mkdir ' + str(folder), shell=True)
if not os.path.exists(str(folder) + "core"):
    subprocess.run('mkdir ' + str(folder) + "core", shell=True)
if not os.path.exists(str(folder) + "flexible"):
    subprocess.run('mkdir ' + str(folder) + "flexible", shell=True)
if not os.path.exists(str(folder) + "unique"):
    subprocess.run('mkdir ' + str(folder) + "unique", shell=True)
core = int(args.l[0])
size = int(args.c[0])
prefix = str(args.p[0])
genomes = str(args.g[0])
min_val = int(args.m[0])
#
#
#
master = []
self = []
sample = []
with open(mmseqs, "r") as infile:
    c = 1 # core
    f = 1 # flexible
    u = 1 # unique
    for name, seq in SFP(infile):
        if seq == "":
            if len(master)/2 >= int(size) and len(list(set(self))) >= int(genomes): # minimum cluster size and minimum unique genomes
                sample_set = list(set(sample))
                sample_len = len(sample_set)
                min_list = []
                min_count = 1
                if min_val > 1:
                    for item in sample_set:
                        min_list.append(sample.count(item))
                    min_count = min(min_list)
                if sample_len >= core and min_count >= min_val: # core proteins
                    i = 0
                    with open(str(folder) + "core/core_cluster_" + str(c) + ".faa", "a") as outfile:
                        while i < len(master):
                            outfile.write(">" + str(master[i]) + "\n" + str(master[i+1]) + "\n")
                            i += 2
                        c += 1
                        master = []
                        self = []
                        sample = []
                elif sample_len == 1: # unique proteins
                    i = 0
                    with open(str(folder) + "unique/unique_cluster_" + str(u) + ".faa", "a") as outfile:
                        while i < len(master):
                            outfile.write(">" + str(master[i]) + "\n" + str(master[i+1]) + "\n")
                            i += 2
                        u += 1
                        master = []
                        self = []
                        sample = []
                elif sample_len > 1 and sample_len < core and min_count >= min_val: # flexible proteins
                    i = 0
                    with open(str(folder) + "flexible/flexible_cluster_" + str(f) + ".faa", "a") as outfile:
                        while i < len(master):
                            outfile.write(">" + str(master[i]) + "\n" + str(master[i+1]) + "\n")
                            i += 2
                        f += 1
                        master = []
                        self = []
                        sample = []
                else: # this should never occur?
                    master = []
                    self = []
                    sample = []
            else:
                master = []
                self = []
                sample = []
        else:
            #if str(name).split(prefix,1)[1] not in master:
            sample.append(str(name).split(prefix,1)[0]) # sample info
            master.append(str(name)) # prot name
            master.append(str(seq)) # prot seq
            self.append(str(name).rsplit("_",1)[0]) # self genome
    #
    #
    # read the last cluster in the file
    if len(master)/2 >= int(size) and len(list(set(self))) >= int(genomes): # minimum cluster size and minimum unique genomes
        sample_set = list(set(sample))
        sample_len = len(sample_set)
        min_list = []
        min_count = 1
        if min_val > 1:
            for item in sample_set:
                min_list.append(sample.count(item))
            min_count = min(min_list)
        if sample_len >= core and min_count >= min_val: # core proteins
            i = 0
            with open(str(folder) + "core/core_cluster_" + str(c) + ".faa", "a") as outfile:
                while i < len(master):
                    outfile.write(">" + str(master[i]) + "\n" + str(master[i+1]) + "\n")
                    i += 2
                c += 1
        elif sample_len == 1: # unique proteins
            i = 0
            with open(str(folder) + "unique/unique_cluster_" + str(u) + ".faa", "a") as outfile:
                while i < len(master):
                    outfile.write(">" + str(master[i]) + "\n" + str(master[i+1]) + "\n")
                    i += 2
                u += 1
        elif sample_len > 1 and sample_len < core and min_count >= min_val: # flexible proteins
            i = 0
            with open(str(folder) + "flexible/flexible_cluster_" + str(f) + ".faa", "a") as outfile:
                while i < len(master):
                    outfile.write(">" + str(master[i]) + "\n" + str(master[i+1]) + "\n")
                    i += 2
                f += 1


#
#
#
