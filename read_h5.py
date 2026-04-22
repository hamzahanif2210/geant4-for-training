import h5py
import numpy as np

file_path = "/scratch/hamza95/train_130k.h5"

with h5py.File(file_path, "r") as f:
    showers = f["showers"]

    # Example: reshape first event
    arr = showers[0]                 # flattened
    reshaped = arr.reshape(-1, 4)    # (N, 4)

    print("Original shape:", arr.shape)
    print("Reshaped:", reshaped.shape)
    print(reshaped[:5])