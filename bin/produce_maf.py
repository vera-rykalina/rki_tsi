#!/usr/bin/env python3

# Import libraries
import pandas as pd
import numpy as np
import sys
import re

infilename = sys.argv[1]
outfilename = sys.argv[2]


# Read .csv file
f = open(infilename, "r")
df = pd.read_csv(f, sep = ",", index_col=False)
f.close()


name1 = infilename.rsplit("/")[-1] # gives a file name.csv
name2 = name1.split("_BaseFreqs_")[0].split("_remap")[0] # gives a sample ID


# Find Max (throughout columns A count: T count)
df["Max"] = df[["A count", "C count", "G count", "T count"]].max(axis=1)

# Find Sum (throughout columns A count: T count)
df["Sum"] = df[["A count", "C count", "G count", "T count"]].sum(axis=1)

# Find Max (throughout columns A count: T count)
df["MAF"] = 1 - df["Max"]/df["Sum"]

# Replace NaNs by zeros 
#df["MAF"] = df["MAF"].replace(np.nan, 0)

# Select only what is needed
df = df.loc[:,["Position in B.FR.83.HXB2_LAI_IIIB_BRU.K03455", "MAF"]]

# Remove raw with "-" at HXB2
df = df[df["Position in B.FR.83.HXB2_LAI_IIIB_BRU.K03455"] != "-"]

# Rename HBX2 column -> position
df.rename({"Position in B.FR.83.HXB2_LAI_IIIB_BRU.K03455": "pos"}, axis=1, inplace=True)

# Rename MAF column -> ID (sample ID)
df.rename({"MAF": name2}, axis=1, inplace=True)

# Prepare a .csv file
df.to_csv(name2 + ".csv", sep=",", index=False, encoding="utf-8")
