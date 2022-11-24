''' cfSimu - a fragment distrubtion simulator '''
import sys
import os
import argparse
from contextlib import redirect_stdout
from io import StringIO

import numpy as np
import deeptools.getFragmentAndReadSize as gfr

from src import Auxillary, parseArgs
from src.estimators import Estimators, Parameters
from test import test_fde


class CfSimu:
    ''' Parent Object Class '''
    def __init__(self, args: argparse.Namespace) -> None:
        ''' Initialize object '''
        self.args = args
        self.help = Auxillary()

        # Output directory
        self.out_dir = self.args.outdir
        self.help.setupDir(self.out_dir)

        # Run name
        self.name = self.args.name
    
        # Setup logging
        self.logger = self.help.logSetup(os.path.join(self.args.outdir, self.name + ".log"))
        self.help.basic_loggings(self.logger, self.args)



class Fde(CfSimu):
    ''' Fragment Distrubtion Estimation '''
    def __init__(self, args: argparse.Namespace) -> None:
        ''' Initialize object '''
        super().__init__(args)
        # Input files
        self.bam_files = self.args.bamfiles
        self.help.checkFileList(self.bam_files)
        # Blacklist file
        self.blacklist_file = self.args.blackListFileName
        if self.blacklist_file:
            self.help.checkFile(self.blacklist_file)
        # Number of processors
        self.nproc = self.args.numberOfProcessors
        # Bin size
        self.bin_size = self.args.binSize
        # Distance between bins
        self.dist_bins = self.args.distanceBetweenBins
        # Max fragment length
        self.maxFragmentLength = self.args.maxFragmentLength

        
    def fde_runner(self) -> None:
        ''' Run Fragment Distrubtion Estimation '''
        # Setup output directory
        self.out_dir = os.path.join(self.out_dir, "fde")
        self.help.setupDir(self.out_dir)
        # Run deeptools
        self.logger.info("Getting fragment size distrubtion")
        self.frag_lens = {}
        for bam in self.bam_files:
            bam_name = os.path.basename(bam.split(".bam")[0])
            # Get fragment size distrubtion from deeptools getFragmentAndReadSize
            bam_pkl = os.path.join(self.out_dir, bam_name + ".frag_len.pkl")
            if os.path.exists(bam_pkl):
                self.logger.info("Found fragment size distrubtion for %s", bam_name)
                self.frag_lens[bam_name] = self.help.read_dict_pkl(bam_pkl)
            else:
                try:
                    with redirect_stdout(StringIO()) as f:
                        frag_len_dict, read_len_dict = gfr.get_read_and_fragment_length(
                                            bam,
                                            return_lengths = True,
                                            blackListFileName = self.args.blackListFileName,
                                            numberOfProcessors = self.args.numberOfProcessors,
                                            verbose = True,
                                            binSize = self.args.binSize,
                                            distanceBetweenBins = self.args.distanceBetweenBins)
                    self.logger.info(f.getvalue().replace('\n', ' '))
                except Exception as e:
                    self.logger.error("Error getting fragment size distrubtion for {}: {}".format(bam_name, e))
                    sys.exit(1)
                self.frag_lens[bam_name] = frag_len_dict
                # frag_len_dict keys:
                # sample_size min mean median max std mad lengths 
                # qtile10 qtile20 qtile25 qtile30 qtile40 qtile60 qtile70 qtile75 qtile80 qtile90 qtile99 
                self.logger.info("Writing fragment size distrubtion for {} to file".format(bam_name))
                self.help.write_dict_pkl(bam_pkl, frag_len_dict)
            
                
    def fde_estimators(self) -> None:
        ''' Compute fragment size distrubtion estimators '''
        self.logger.info("Computing fragment size distrubtion estimators")
        for name, fragment_dict in self.frag_lens.items():
            self.logger.info("Computing fragment size distrubtion estimators for {}".format(name))
            # Estimate
            fde_params = Parameters(name, self.logger, self.out_dir, self.maxFragmentLength, fragment_dict)
            estimator = Estimators(fde_params)
            
def run():
    '''Run cfSimu'''
    # Parse command line arguments
    args = parseArgs()
    # Switch to subcommand
    subcmd = args.subparser_name
    if subcmd == "fde":
        # Run Fragment Distrubtion Estimation
        fde = Fde(args)
        fde.fde_runner()
        fde.fde_estimators()
    elif subcmd == "tester":
        # Run test
        test_fde()
    else:
        # Run cfSimu
        simu = CfSimu(args)

    
if __name__ == '__main__':
    run()
    