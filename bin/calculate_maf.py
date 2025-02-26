#!/usr/bin/env python3

# Import libraries
import argparse
import pandas as pd
import numpy as np


# Global variables 
_args = None

def initialise():
    '''
    Parse command-line arguments.
    '''
    global _args
    parser = argparse.ArgumentParser( description="Calculation of minor allele frequency (MAF) as MAF = (1-(A+C+G+T))/(A+C+G+T)" )
    parser.add_argument( "-b", "--basefreqspath", required=True, help="A comma-separated file from shiver: BaseFreqWithHXB2." )
    parser.add_argument( "-r", "--hxb2", required=True, help="A comma-separated HXB2 reference file" )
    parser.add_argument( "-o", "--outpath", required=True, help="Filename for saving calculated MAF." )
    _args = parser.parse_args()
    return

def load_shiver_basefreqs(fpath):
    ''' Load shiver output to calculate cumumative minor allele frequencies. '''
    basefreqs = pd.read_csv(fpath, sep=",", index_col=False)
    print("Loaded shviver BaseFreqs with HXB2 data, shape={}".format(basefreqs.shape))
    # Extract sample ID from a shiver input file
    sample_id = fpath.rsplit("/")[-1].split("_BaseFreqs_")[0].split("_remap")[0]
    return basefreqs, sample_id

def load_reference(fpath):
    ''' Load HXB2 reference. '''
    hxb2 = pd.read_csv(fpath, sep=",", index_col=False)
    print("Loaded shviver BaseFreqs with HXB2 data, shape={}".format(hxb2.shape))
    return hxb2

def calculate_maf(basefreqs, sample_id ):
    ''' Calculate MAF. '''
    maf = basefreqs.copy()
    # Find Max (throughout columns A count: T count)
    maf["Max"] = maf[["A count", "C count", "G count", "T count"]].max(axis=1)
    # Find Sum (throughout columns A count: T count)
    maf["Sum"] = maf[["A count", "C count", "G count", "T count"]].sum(axis=1)
    # Find Max (throughout columns A count: T count)
    maf["MAF"] = 1 - maf["Max"]/maf["Sum"]
    # Remove raws with "-" at HXB2
    maf = maf[maf["Position in B.FR.83.HXB2_LAI_IIIB_BRU.K03455"] != "-"]
    # Select only what is needed
    maf = maf.loc[:,["Position in B.FR.83.HXB2_LAI_IIIB_BRU.K03455", "MAF"]]
    # Rename HBX2 column -> pos and MAF column -> sample ID (variable)
    maf.rename({"Position in B.FR.83.HXB2_LAI_IIIB_BRU.K03455": "pos", "MAF": sample_id}, axis=1, inplace=True)
    maf["pos"] = maf["pos"].astype(np.int64)
    return maf

def transform_maf(hxb2, maf):
    ''' Reshape MAF file '''
    ref_maf = hxb2.merge(maf, on="pos", how="left")
    # Remove a sanity check column
    ref_maf.drop(["HXB2 base"], axis=1, inplace=True)
    # Set index to "pos" columns
    ref_maf = ref_maf.set_index("pos")
    # Transpose df
    transposed_maf = ref_maf.T
    # Reset index
    transposed_maf = transposed_maf.reset_index()
    return transposed_maf

def save_maf(transposed_maf, outf):
    ''' Write output. '''
    transposed_maf.to_csv(outf, sep=",", header = True, index=False, encoding="utf-8")
    print("Output saved as {}.".format(outf))
    print(transposed_maf.head())

def main():
    ''' Run calculations. '''
    basefreqs, sample_id = load_shiver_basefreqs(_args.basefreqspath)
    hxb2 = load_reference(_args.hxb2)
    maf = calculate_maf(basefreqs, sample_id)
    transposed_maf = transform_maf(hxb2, maf)
    save_maf(transposed_maf, _args.outpath)

if __name__ == '__main__':
    initialise()
    main()