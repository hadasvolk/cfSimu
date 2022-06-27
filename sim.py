# @plac.pos('genome', "FASTA reference sequence")
# @plac.pos('read1', "FASTQ file for first in pair (read1.fq)")
# @plac.pos('read2', "FASTQ file for second in pair (read2.fq)")
# @plac.opt('err', "the base error rate", abbrev='e', type=float)
# @plac.opt('mut', "rate of mutations", abbrev='r', type=float)
# @plac.opt('frac', "fraction of indels", abbrev='R', type=float)
# @plac.opt('ext', "probability an indel is extended", abbrev='X', type=float)
# @plac.opt('dist', "outer distance between the two ends", abbrev='D', type=int)
# @plac.opt('stdev', "standard deviation", abbrev='s', type=int)
# @plac.opt('L1', "length of the first read", abbrev='1', type=int)
# @plac.opt('L2', "length of the second read", abbrev='2', type=int)
# @plac.opt('num', "number of read pairs", abbrev='N', type=int)
# @plac.opt('seed', "seed for the random generator", abbrev='S', type=int)
# @plac.opt('amb', "disregard if the fraction of ambiguous bases higher than FLOAT", abbrev='A', type=float)
# @plac.flg('fixed', "each chromosome gets N sequences", abbrev='f')
# @plac.flg('version', "print version number", abbrev='v')

import sys, os, logging
import yaml
import random
import pickle
import argparse
import pandas as pd
from pywgsim import wgsim
from Bio import SeqIO, bgzf
from ray.util.multiprocessing import Pool

REF = '/storage/data/resources/reference_genome/hg38/canonical/hg38.canonical.fa'


def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', '-n', help='name of the simulation', type=str, required=True)
    parser.add_argument('--ref', '-r', help='Genome FASTA file', type=str, default=REF, required=False)
    parser.add_argument('--lengths', '-l', help='Lengths file', type=str, default='maternal_lengths.csv', required=False)
    parser.add_argument('--n_seqs', '-ns', help='Number of sequences', type=int, default=148576028, required=False)
    parser.add_argument('--n_workers', '-nw', help='Number of workers', type=int, default=1, required=False)
    parser.add_argument('--outdir', '-o', help='Output directory', type=str, default='simulations', required=False)
    parser.add_argument('--fastq_list', '-fq', help='List of fastq files', type=str, default='fastq_list.txt', required=False)
    return parser.parse_args()


class ReadGenerator:
    def __init__(self, args) -> list:
        self.name = args.name
        self.lengths_file = args.lengths
        self.n_seqs = args.n_seqs
        self.ref = args.ref
        self.n_workers = args.n_workers
        self.MAX_READS = 100000
        
        self.out = args.outdir
        self.fastq_list = args.fastq_list

        try:
            os.makedirs(self.out, exist_ok=True)
        except OSError:
            print("Error: Creating output directory")
            exit()

        self.gendf()

        all_records = []
        all_records_pkl = os.path.join(self.out, '{}.all_records.pkl'.format(self.name))
        if os.path.exists(all_records_pkl):
            with open(all_records_pkl, 'rb') as f:
                all_records = pickle.load(f)
        else:
            for l in self.lengths_df['length'].unique():
                i = 0
                mod = divmod(self.lengths_df[self.lengths_df['length'] == l]['n_reads'].values[0], self.MAX_READS)
                if mod[0] != 0:
                    for i in range(mod[0]):
                        all_records.append([l, self.MAX_READS, i])
                    all_records.append([l, mod[1], i+1])
                else:
                    all_records.append([l, mod[1], 0])
            pickle.dump(all_records, open(all_records_pkl, 'wb'))
        print("number of records:",len(all_records))
        str_all_rec = [','.join(map(str, l)) for l in all_records]
        with Pool(processes = self.n_workers) as pool:
            pool.map(self.exec, str_all_rec)
        # for l in str_all_rec:
        #     self.exec(l)

        a = {}
        for record in all_records:
            r1 = "{}/read1.{}.{}.fq.gz".format(self.out, record[0], record[2])
            r2 = "{}/read2.{}.{}.fq.gz".format(self.out, record[0], record[2])
            a["{}-{}".format(record[0], record[2])] = [r1, r2] 
        print("number of records:",len(all_records))
        with open(self.fastq_list, 'wb') as handle:
            pickle.dump(all_records, handle)


    def gendf(self) -> None:
        df = pd.read_csv(self.lengths_file)
        df = df[(df['length'] > 49) & (df['length'] < 501)]
        df['norm'] = df['COUNT(length)'] / df['COUNT(length)'].sum()
        df['n_reads'] = df['norm'] * self.n_seqs
        df['n_reads'] = df['n_reads'].astype(int)
        self.lengths_df = df


    def exec(self, record):
        print(record)
        dist = int(record.split(',')[0])
        N = int(record.split(',')[1])
        s = random.randint(1, 500)
        mutation_rate = 0
        print('dist:\t{}\t\tN:\t{}\t\tseed:\t{}'.format(dist, N, s))
        r1 = "{}/read1.{}.{}.fq".format(self.out, dist, record.split(',')[2])
        r2 = "{}/read2.{}.{}.fq".format(self.out, dist, record.split(',')[2])
        if dist < 151:
            L = dist
        else:
            L = 151
        
        wgsim.core(r1=r1, r2=r2, 
                    ref=self.ref, 
                    err_rate=0.006, 
                    mut_rate=mutation_rate, 
                    indel_frac=0.15, 
                    indel_ext=0.25, 
                    max_n=0.05, 
                    is_hap=0, 
                    N=N,  
                    dist=dist, 
                    stdev=50, 
                    size_l=L, 
                    size_r=L, 
                    is_fixed=0, 
                    seed=s)
        
        seqs1 = []
        seqs2 = []        
        def rename(file, is_read1 = True):
            if os.stat(file).st_size == 0: return
            with open(file, 'r') as f:
                for seq in SeqIO.parse(f, "fastq"):
                    seq.name = None
                    if is_read1:
                        seq.id = "{}|{}|{}/{}".format(dist, record.split(',')[2], seq.id, 1)
                        seqs1.append(seq)
                    else:
                        seq.id = "{}|{}|{}/{}".format(dist, record.split(',')[2], seq.id, 2)
                        seqs2.append(seq)

        rename(r1) 
        rename(r2, is_read1 = False)
        try:
            os.remove(r1)
            os.remove(r2)
        except OSError:
            print("Error: removing files {} and {}".format(r1, r2))

        for i in zip([r1, r2], [seqs1, seqs2]): self.write_fastq_compress(i)

    
    def write_fastq_compress(self, seqs):
        with bgzf.BgzfWriter("{}.gz".format(seqs[0]), "wb") as outgz:
            SeqIO.write(sequences=seqs[1], handle=outgz, format="fastq")


if __name__ == '__main__':   
    ReadGenerator(parseArgs())