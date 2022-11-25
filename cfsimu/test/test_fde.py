import sys
import os
import pathlib
import argparse

import cfsimu

def test_fde():
    ''' Test Fragment Distrubtion Estimation '''
    args = argparse.Namespace(**{
        "subparser_name": "fde",
        "bamfiles": [os.path.realpath(os.path.join("data", "bams", "S36_1.sorted.mdup.bam"))],
        "name": "test",
        "outdir": "cfsimu/test/output",
        "numberOfProcessors": 40,
        "binSize": 10,
        "distanceBetweenBins": 10000,
        "blackListFileName": None,
        "maxFragmentLength": 1000,
        "loss": "linear",
        "f_scale": 1.0
    })

    fde = cfsimu.Fde(args)
    fde.fde_runner()
    fde.fde_estimators()
