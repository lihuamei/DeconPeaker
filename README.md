DeconPeaker
===================================================

`DeconPeaker`: a deconvolution method to estimate cell type proportions in chromatin accessibility data (ATAC-Seq), as well as gene expression data (RNA-Seq & Microarray).
![DeconPeaker\_pipeline](pipeline.png)

How to use `DeconPeaker`?
---------------------
DeconPeaker's code is a mix of Python3.8 and R(4.0), which requires the following dependencies.
* Python3.8:
	* Numpy
	* Scipy
	* Pandas
	* bx
	* Matplotlib
	* rpy2
* R4.0:
	* pls
	* transport
	* colorRamps
	* MASS
* Other tools (when excute preprocess and simulation steps):
	* bedtools
	* samtools
	* featureCounts

Installation of DeconPeaker
--------------------------
```bash
git clone git@github.com:lihuamei/DeconPeaker.git
cd DeconPeaker
pip install .

cd ..
```

```python
deconPeaker --help
usage: deconPeaker [-h] {preprocess,findctsps,deconvolution,simulation} ...

deconPeaker - a deconvolution model to identify cell types based on chromatin accessibility in ATAC-Seq data of mixture samples.

positional arguments:
  {preprocess,findctsps,deconvolution,simulation}
    preprocess          Create chromatin accessibility profile for pure samples, this step is the basis for subsequent specific cell type identification and mixture deconvolution. Note: This step only support Linux system, and if
                        you have a large sample size for pure cells, please ensure enough sufficient memory and hard storage space for program to run normally.
    findctsps           Find cell type specific peaks/genes accross pure samples, different pure cell samples require replicates as input.
    deconvolution       Based on pure cell profile information, robust regression deconvolution strategy was used to estimate the proportion of possible cell types in the mixed samples.
    simulation          Simulate mixed samples with different proportions of cells. [Note] The method is proportional random sampling of reads from different cell types of BAM/BED files.

optional arguments:
  -h, --help            show this help message and exit
INFO  @ Tue, 23 Apr 2024 15:09:32: Embedded R ended.
INFO  @ Tue, 23 Apr 2024 15:09:32: Embedded R already ended.

```

Show an exmaple

```bash
deconPeaker deconvolution --lib-strategy=ATAC-Seq --mixture=DeconPeaker/test/examples/ATAC-Seq/GSE74912_Corces_MR_synthetic_mixture_counts_data.xls --pure=DeconPeaker/test/examples/ATAC-Seq/GSE74912_Corces_MR_pure_readcounts_signature_matrix.xls --format=TABLE --pvalue=FALSE --outdir=DeconPeaker/results/GSE74912_Corces_MR

INFO  @ Tue, 23 Apr 2024 14:49:47: Deconvolving......
[1] "[WARN] Lambda = 0, Log transfer for SIMPLS"
[1] "[WARN] Lambda = 0, Log transfer for SIMPLS"
[1] "[WARN] Lambda = 0, Log transfer for SIMPLS"
INFO  @ Tue, 23 Apr 2024 14:50:17: Showing deconPeaker results:
                 Ery       CLP      CD8T       GMP        NK       CMP      CD4T         B       MEP       MPP      LMPP       HSC      MONO  Rsquared     RMSEP  P.value
Sample_171  0.069005  0.002305  0.052112  0.067137  0.262002  0.013519  0.286919  0.080115  0.049956  0.000000  0.047004  0.065608  0.004319  0.866373  0.365072   9999.0
Sample_177  0.018484  0.111797  0.000000  0.047688  0.263781  0.066290  0.007585  0.049847  0.056509  0.041708  0.218438  0.074697  0.043177  0.928823  0.266408   9999.0
Sample_159  0.000289  0.286506  0.000000  0.090313  0.032724  0.016271  0.245973  0.071600  0.000000  0.191758  0.014806  0.046406  0.003353  0.964885  0.187197   9999.0
Sample_139  0.006446  0.008090  0.000000  0.000000  0.000000  0.011153  0.909165  0.000000  0.006506  0.000000  0.022982  0.029292  0.006367  0.898426  0.318128   9999.0
Sample_86   0.382248  0.165546  0.021309  0.015785  0.003639  0.090452  0.240281  0.024951  0.025254  0.022799  0.000000  0.007111  0.000624  0.952080  0.218717   9999.0
...              ...       ...       ...       ...       ...       ...       ...       ...       ...       ...       ...       ...       ...       ...       ...      ...
Sample_134  0.002720  0.000000  0.019477  0.019017  0.000000  0.017883  0.000000  0.000000  0.000000  0.082663  0.012248  0.834613  0.011378  0.901643  0.313398   9999.0
Sample_1    0.002720  0.000000  0.019477  0.019017  0.000000  0.017883  0.000000  0.000000  0.000000  0.082663  0.012248  0.834613  0.011378  0.901719  0.313285   9999.0
Sample_84   0.019543  0.005681  0.000000  0.008541  0.233237  0.324421  0.211760  0.000000  0.030085  0.000000  0.124828  0.041904  0.000000  0.911399  0.297173   9999.0
Sample_175  0.000891  0.059974  0.021029  0.000000  0.196285  0.253763  0.208002  0.000000  0.004527  0.064388  0.047768  0.136915  0.006459  0.908695  0.301642   9999.0
Sample_29   0.000000  0.044709  0.002532  0.141323  0.000000  0.028472  0.000000  0.000000  0.000000  0.000000  0.674772  0.064291  0.043901  0.980558  0.139116   9999.0

[195 rows x 16 columns]
INFO  @ Tue, 23 Apr 2024 14:50:17: Elapsed time is 29.629956245422363 seconds
INFO  @ Tue, 23 Apr 2024 14:50:17: Embedded R ended.
INFO  @ Tue, 23 Apr 2024 14:50:17: Embedded R already ended.
```
-------------------
It is necessary to mention that there are four parts in DeconPeaker, and details as follows:
* preprocess: only supports on Linux system
* findctsps: supports on Windows and Linux systems
* deconvolution: supports on Windows and Linux systems
* simulation: only supports on Linux system

More Information
--------------------
Please see [Tutorial](https://lihuamei.github.io//DeconPeaker/test/DeconPeak_demo.html).

Citation
---------------------
Please cite the publication: ***Li H, Sharma A, Luo K, et al. DeconPeaker, a deconvolution model to identify cell types based on chromatin accessibility in ATAC-Seq data of mixture samples[J]. Frontiers in genetics, 2020, 11: 392.***<br>
