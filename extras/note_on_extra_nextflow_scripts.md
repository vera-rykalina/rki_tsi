## Notes

This folder contains two additional scripts:

- **hivtime_singles.nf**  
  A convenient script for running individual samples to enrich the training dataset. It organizes the result folder, ensuring that file names are tagged with the sample ID for easy identification.

- **fastq_pol_filter.nf**  
  This script extracts the **pol** region from full-length FASTQ files, generating separate **pol** FASTQ files.

## NB!  
If you want to use any of these scripts, please copy them to the main pipeline folder (one level up). They are designed to run from the **rki_tsi/** folder, and will not work if run from elsewhere.
