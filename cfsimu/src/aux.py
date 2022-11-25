import argparse
import sys 
import os
import logging
import pickle as pkl

import pandas as pd
import numpy as np


class Auxillary:
    ''' Helper functions '''
    
    @staticmethod
    def setupDir(dir: str) -> None:
        ''' Setup if directory exists '''
        try:
            os.makedirs(dir, exist_ok=True)
        except OSError as error:
            logging.error("Error: unable to create directory: {}\n".format(dir, error))
            sys.exit(1)

            
    @staticmethod
    def checkFile(file: str) -> None:
        ''' Check if file exists '''
        if not os.path.isfile(file):
            sys.exit("File {} does not exist".format(file))

            
    @staticmethod
    def checkFileList(file_list: list) -> None:
        ''' Check if list of files exists '''
        for file in file_list:
            if not os.path.isfile(file):
                logging.error("File {} does not exist".format(file))
                sys.exit(1)
    
                
    def logSetup(self, logfile: str) -> logging.Logger:
        ''' Setup logging '''
        logging.basicConfig(level=logging.INFO, 
                            format="%(asctime)s %(message)s", 
                            datefmt="%d/%m/%Y %H:%M:%S",
                            handlers=[logging.FileHandler(logfile), logging.StreamHandler()])
        return logging.getLogger('cfSimu')
    

    def basic_loggings(self, logger: logging.Logger, args: argparse.Namespace) -> None:
        ''' Log versions '''
        logger.info("Starting cfSimu - {}".format(args.subparser_name))
        logger.info("Working directory: {}".format(os.getcwd()))
        for key, value in vars(args).items():
            logger.info("Argument: {} = {}".format(key, value))
        # logger.info("Python version: {}".format(sys.version.split()[0]))
        # logger.info("Pandas version: {}".format(pd.__version__))
        # logger.info("Numpy version: {}".format(np.__version__))
    
    
    def write_dict_pkl(self, out: str, dict: dict) -> None:
        ''' Write dict to pickle '''
        try:
            with open(out, 'wb') as handle:
                pkl.dump(dict, handle, protocol=pkl.HIGHEST_PROTOCOL)
        except OSError as error:
            logging.error("Error: unable to write to file: {}\n".format(out, error))
            sys.exit(1)
    
            
    def read_dict_pkl(self, file: str) -> dict:
        ''' Read dict from pickle '''
        try:
            with open(file, 'rb') as handle:
                dict = pkl.load(handle)
        except OSError as error:
            logging.error("Error: unable to read file: {}\n".format(file, error))
            sys.exit(1)
        return dict
    

    def write_df_tsv(self, out: str, df: pd.DataFrame) -> None:
        ''' Write dataframe to tsv '''
        try:
            df.to_csv(out, sep='\t', index=False)
        except OSError as error:
            logging.error("Error: unable to write to file: {}\n".format(out, error))
            sys.exit(1)
    

        
def parseArgs() -> argparse.Namespace:
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser(prog="cfSimu")
    subparsers = parser.add_subparsers(help="sub-command help", 
                            dest="subparser_name")

    # create the parser for the Fragment Distrubtion Estimation command
    parser_fde = subparsers.add_parser("fde", 
                            help="Fragment Distrubtion Estimation")
    
    parser_fde.add_argument("--bamfiles", '-b',
                            help="List of BAM files to process",
                            nargs='+',
                            required=True)
    
    parser_fde.add_argument("--name", '-n', 
                            help="name of the simulation", 
                            type=str,
                            required=True)

    parser_fde.add_argument("--outdir", '-o', 
                            help="Output directory (Default: %(default)s)", 
                            type=str, 
                            default="cfSimu_out", 
                            required=False)

    # DeepTools arguments
    parser_fde.add_argument("--numberOfProcessors", '-p',
                            help="Number of processors to use. The default is "
                            "to use 1. (Default: %(default)s)",
                            metavar="INT",
                            type=int,
                            default=1,
                            required=False)
    
    parser_fde.add_argument("--binSize", "-bs",
                            metavar="INT",
                            help="Length in bases of the window used to sample the genome. (Default: %(default)s)",
                            default=1000,
                            type=int)
    
    parser_fde.add_argument("--distanceBetweenBins", "-dn",
                            metavar="INT",
                            help="To reduce the computation time, not every possible genomic "
                            "bin is sampled. This option allows you to set the distance "
                            "between bins actually sampled from. Larger numbers are sufficient "
                            "for high coverage samples, while smaller values are useful for "
                            "lower coverage samples. Note that if you specify a value that "
                            "results in too few (<1000) reads sampled, the value will be "
                            "decreased. (Default: %(default)s)",
                            default=1000000,
                            type=int)
    
    parser_fde.add_argument("--blackListFileName", "-bl",
                            help="A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered.",
                            metavar="BED FILE",
                            required=False)
    
    parser_fde.add_argument("--maxFragmentLength", "-ml",
                            metavar="INT",
                            help="Maximum fragment length to consider. (Default: %(default)s)",
                            default=1000,
                            type=int)

    # Least squares arguments
    parser_fde.add_argument("--loss", "-l",
                            help="Loss function to use in the least_square method. (Default: %(default)s)",
                            choices=["soft_l1", "huber", "cauchy", "arctan"],
                            default="linear",
                            required=False)

    parser_fde.add_argument("--f_scale", "-fs",
                            help="Scale factor to multiply the loss function with. (Default: %(default)s)",
                            type=float,
                            default=1.0,
                            required=False)
    
    

    # Create the parser for the tester module
    parser_tester = subparsers.add_parser("tester",
                            help="Tester module")
    
    # parser.add_argument('--ref', '-r', help='Genome FASTA file', type=str, default=REF, required=False)
    # parser.add_argument('--lengths', '-l', help='Lengths file', type=str, required=True)
    # parser.add_argument('--n_seqs', '-ns', help='Number of sequences', type=int, default=148576028, required=False)
    # parser.add_argument('--n_workers', '-nw', help='Number of workers', type=int, default=1, required=False)
    # parser.add_argument('--outdir', '-o', help='Output directory', type=str, default='simulations', required=False)
    # parser.add_argument('--fastq_dict', '-fq', help='Pickle dict of fastq files produced', type=str, default='fastq_dict', required=False)

    if len(sys.argv)==1:
        sys.stderr.write("\nNo subcommand given\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()
    
