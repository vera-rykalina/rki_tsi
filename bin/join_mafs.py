#!/usr/bin/env python3

# Import libraries
import pandas as pd
import sys


# Create a list of dfs and assign HXB2_refdata.csv as a separate df
dfs = []
for infilename in sys.argv[1:]:
    name = infilename.rsplit("/")[-1] # gives a file name name.csv
    list_of_substrings = name.split("_")
    if "HXB2" in list_of_substrings:
        df = pd.read_csv(infilename, sep = ",", index_col=False)
    else:
        dfs.append(pd.read_csv(infilename, sep = ",", index_col=False))


# Join each df in dfs to the HXB2_refdata.csv df
for df_ in dfs[0:]:
    df = df.merge(df_, on="pos", how="left")

# Remove a sanity check column
df.drop(["HXB2 base"], axis=1, inplace=True)

# Fill empty cells with zeros (a need of this step should be double ckecked!!!)
#df.fillna(0, inplace = True,  downcast="infer")      

# Set index to "pos" columns
df = df.set_index("pos")

# Transpose df
df = df.T

# Reset index
df = df.reset_index()
print(df)

# Create a clean csv file
df.to_csv("joined" + "_MAF" + ".csv", sep=",", header = True, index = False, encoding="utf-8")
