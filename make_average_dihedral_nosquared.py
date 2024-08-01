import numpy as np


def load_data(filename):
    """
    Load data from a text file. Assumes that the file has columns of numeric data.
    """
    try:
        data = np.loadtxt(filename, delimiter=',', usecols=(1, 2))
        return data
    except Exception as e:
        print(f"Error loading data from file: {e}")
        return None

def compute_statistics(data):
    """
    Compute the mean and standard deviation for each column in the data.
    """
    if data is None:
        return None
    
    means = np.mean(data, axis=0)
    std_devs = np.std(data, axis=0)
    
    return means, std_devs

def main(filename):
    data = load_data(filename)
    
    if data is not None:
        means, std_devs = compute_statistics(data)
        
        for i, (mean, std_dev) in enumerate(zip(means, std_devs), start=1):
            print(f"Column {i}: Mean = {mean}, Standard Deviation = {std_dev}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
    else:
        main(sys.argv[1])

