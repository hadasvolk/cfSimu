import logging

import numpy as np
import pandas as pd
import scipy as sp
from scipy.signal import argrelextrema
from scipy.optimize import least_squares
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

    # Loss function for least squares
    loss: str = "linear"
    # Scale factor for the loss function
    f_scale: float = 0.1

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
        return self.maxima


class Estimators:
    def __init__(self, params: Parameters) -> None:
        ''' Initialize estimators attributes '''
        self.params = params
        self.params.logger.info("")
        self.params.logger.info("Initailized Estimators object for {}".format(params.name))
    
        
    def exec_leastsq(self, func: callable, x0: np.array, args: tuple) -> sp.optimize.optimize.OptimizeResult:
        ''' Execute least squares '''
        loss = self.params.loss
        f_scale = self.params.f_scale
        self.params.logger.info("Executing least squares")
        self.params.logger.info("Estimate function: {}".format(func.__name__))
        self.params.logger.info("Initial guess: {}".format(x0))
        self.params.logger.info("Loss function: {}".format(loss))
        self.params.logger.info("Scaling factor: {}".format(f_scale))
        # Execute least squares
        res = least_squares(func, 
                            x0, 
                            verbose = 0, 
                            loss = loss, 
                            f_scale = f_scale,
                            bounds = (0, np.inf), 
                            args=args)
        # Print to log
        self.params.logger.info("Least squares message: {}".format(res.message))
        self.params.logger.info("Least squares cost: {}".format(res.cost))
        self.params.logger.info("Least squares success: {}".format(res.success))
        self.params.logger.info("Least squares optimization status: {}".format(res.status))
        self.params.logger.info("Least squares optimization reason: {}".format(res.optimality))
        self.params.logger.info("Least squares optimization number of iterations: {}".format(res.nfev))
        self.params.logger.info("Least squares optimization number of function evaluations: {}".format(res.njev))
        self.params.logger.info("Least squares optimization number of jacobian evaluations: {}".format(res.njev))
        self.params.logger.info("Least squares optimization x: {}".format(list(res.x)))
        
        # Return result
        return res

    
    def double_gaussian_residual(self, x0: np.array, y: pd.Series, t: pd.Series) -> np.array:
        ''' Double gaussian residual '''
        # Compute residual
        return self.double_gaussian(x0, t) - y
    

    def double_gaussian(self, x0: np.array, t: pd.Series) -> np.array:
        ''' Double gaussian '''
        # Unpack x
        A1, mu1, sigma1, A2, mu2, sigma2 = x0
        return  A1 * np.exp( - (t - mu1) ** 2.0 / (2.0 * sigma1 ** 2.0)) \
                + A2 * np.exp( - (t - mu2) ** 2.0 / (2.0 * sigma2 ** 2.0))

    
    def exec_double_gaussian(self) -> sp.optimize.optimize.OptimizeResult:
        ''' Double Gaussian Estimation '''
        # Sigma values for double gaussian
        SIGMA_1 = 10.0
        SIGMA_2 = 20.0

        self.params.logger.info("")
        self.params.logger.info("Running double Gaussian estimation for {}".format(self.params.name))
        # Fragment lengths
        self.x = self.params.sub_frag_lens_df["frag_len"]
        # Fragment length counts
        self.y = self.params.sub_frag_lens_df["count"]
        # Get max peaks
        max_peaks = self.params.get_max_peaks()
        peak_1_idx = sorted(max_peaks.keys())[0]
        peak_2_idx = sorted(max_peaks.keys())[1]
        self.params.logger.info("Peak 1 index: {}, counts {}".format(peak_1_idx, max_peaks[peak_1_idx]))
        self.params.logger.info("Peak 2 index: {}, counts {}".format(peak_2_idx, max_peaks[peak_2_idx]))

        # Least squares fit. Starting values found by inspection.
        return self.exec_leastsq(self.double_gaussian_residual, 
                    [max_peaks[peak_1_idx], peak_1_idx, SIGMA_1, max_peaks[peak_2_idx], peak_2_idx, SIGMA_2],
                    (self.y, self.x))
        
    
    def plot_estimations(self, fits: dict) -> None:
        ''' Plot estimations '''
        self.params.logger.info("")
        self.params.logger.info("Plotting estimations for {}".format(self.params.name))

        # Plot
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))
        for name, fit in fits.items():
            covar = [float(c) for c in list(fit.x)]
            ax.plot(self.x, self.double_gaussian(covar, self.x), label=name, c='r')

        ax.scatter(self.x, self.y, s=1)
        ax.set_xlabel("Fragment Length")
        ax.set_ylabel("Count")
        ax.set_title("Fragment Length Distribution")
        plt.show()

