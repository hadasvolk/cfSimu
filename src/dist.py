import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn

if __name__ == '__main__':
    # Read in the data
    data = pd.read_csv('data/creditcard.csv')
    # Plot histograms of each parameter
    data.hist(figsize = (20, 20))
    plt.show()
    # Determine number of fraud cases in dataset
    fraud = data[data['Class'] == 1]
    valid = data[data['Class'] == 0]