#!/usr/bin/env python3

# Import libraries
import argparse
import pandas as pd


# Global variables 
_args = None

def initialise():
    '''
    Parse command-line arguments.
    '''
    global _args
    parser = argparse.ArgumentParser( description="Calculation of minor allele frequency (MAF) as MAF = (1-(A+C+G+T))/(A+C+G+T)" )
    parser.add_argument( "-b", "--basefreqspath", required=True, help="A comma-separated file from shiver: BaseFreqWithHXB2." )
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
    return maf

def save_maf(maf, outf):
    ''' Write output. '''
    maf.to_csv(outf, sep=",", index=False, encoding="utf-8")
    print("Output saved as {}.".format(outf))

def main():
    ''' Run calculations. '''
    basefreqs, sample_id = load_shiver_basefreqs(_args.basefreqspath)
    maf = calculate_maf(basefreqs, sample_id)
    save_maf(maf, _args.outpath)

if __name__ == '__main__':
    initialise()
    main()