import logging

import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from dataclasses import dataclass

from src.persistence1d import RunPersistence


@dataclass
class Parameters:
    ''' Parameters for Estimators '''
    # Sample name
    name: str
    # Logger
    logger: logging.Logger
    # Output directory
    out_dir: str
    # Maximum fragment length
    maxFragmentLength: int
    # Fragment lengths dictonary
    frag_lengths_dict: dict
    # Persistence threshold
    presistence_threshold: int = 10000

    
    def __post_init__(self) -> None:
        # Fragment length array
        self.frag_lengths = self.frag_lengths_dict["lengths"]
        # Sample size
        self.sample_size = self.frag_lengths_dict["sample_size"]
        # Fragment length mean
        self.frag_length_mean = self.frag_lengths_dict["mean"]
        # Fragment length standard deviation
        self.frag_length_std = self.frag_lengths_dict["std"]
        # Fragment length median
        self.frag_length_median = self.frag_lengths_dict["median"]
        # Fragment length mad
        self.frag_length_mad = self.frag_lengths_dict["mad"]
        # Fragment length min
        self.frag_length_min = self.frag_lengths_dict["min"]
        # Fragment length max
        self.frag_length_max = self.frag_lengths_dict["max"]

        # Compute fragment length distrubtion in pandas dataframe and trimmed for max fragment length
        unique, counts = np.unique(self.frag_lengths, return_counts=True)
        self.frag_lengths_df = pd.DataFrame({"frag_len": unique, "count": counts}, dtype=int).sort_values(by="frag_len")
        self.sub_frag_lens_df = self.frag_lengths_df[
                (self.frag_lengths_df["frag_len"] > 0) & (self.frag_lengths_df["frag_len"] < self.maxFragmentLength)]
        # Zero fragment length count
        self.zero_frag_len_count = self.frag_lengths_df["count"].values[0]
    
    
    def get_max_peaks(self) -> dict:
        ''' Get n max peaks '''
        self.logger.info("Getting max peaks using persistence1d")
        # Max frgament occurance
        # This simple call ompute the extrema of the given data and their persistence.
        ExtremaAndPersistence = RunPersistence(self.sub_frag_lens_df["count"].to_numpy())
        # Keep only those extrema with a persistence larger than the threshold.
        Filtered = [t for t in ExtremaAndPersistence if t[1] > self.presistence_threshold]
        # Sort the list of extrema by persistence.
        Sorted = sorted(Filtered, key=lambda ExtremumAndPersistence: ExtremumAndPersistence[0])
        # Get the extrema with the highest persistence.
        self.maxima = {}
        # Print to log
        for i, E in enumerate(Sorted):
            if i % 2 != 0:
                self.maxima[E[0]] = self.sub_frag_lens_df["count"].to_numpy()[E[0]]
                self.logger.info("Maximum at index {} with persistence {} and data value {}".format( 
                                                        E[0],
                                                        E[1], 
                                                        self.sub_frag_lens_df["count"].to_numpy()[E[0]]))


class Estimators:
    def __init__(self, params: Parameters) -> None:
        ''' Initialize estimators attributes '''
        params.logger.info("Initailized Estimators object for {}".format(params.name))

    
    def double_gauss(self) -> None:
        ''' Double Gaussian Estimation '''
        self.logger.info("Running double Gaussian estimation for {}".format(self.name))
        # Max frgament occurance
        #~ This simple call ompute the extrema of the given data and their persistence.
        ExtremaAndPersistence = RunPersistence(self.sub_frag_lens_df["count"].to_numpy())

        #~ Keep only those extrema with a persistence larger than 10000.
        Filtered = [t for t in ExtremaAndPersistence if t[1] > 10000]

        #~ Sort the list of extrema by persistence.
        Sorted = sorted(Filtered, key=lambda ExtremumAndPersistence: ExtremumAndPersistence[0])

        #~ Print to log
        for i, E in enumerate(Sorted):
            if i % 2 != 0:
                self.logger.info("Maximum at index {} with persistence {} and data value {}".format( 
                                                        E[0],
                                                        E[1], 
                                                        self.sub_frag_lens_df["count"].to_numpy()[E[0]]))
            
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        ax.scatter(self.sub_frag_lens_df["frag_len"], self.sub_frag_lens_df["count"], s=1)
        ax.set_xlabel("Fragment Length")
        ax.set_ylabel("Count")
        ax.set_title("Fragment Length Distribution")
        plt.show()

