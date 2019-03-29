# This is the preprocessing step of DeconPeaker, containing two processes;
# 1) generate a list non-redundant peaks with high confidence;
# 2) count the fragments across all non-redundant peaks based BAM of aligned files.

python deconPeaker.py preprocess --pure=test/GSE74912_ATAC_pure_cells.yaml --thread=10 --prefix=Test --outdir=results


