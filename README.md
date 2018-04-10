# Introduction  
After extensively comparing different WGS-based CNA tools, including [WISECONDOR](https://github.com/VUmcCGP/wisecondor),
[QDNAseq](https://github.com/ccagc/QDNAseq), [CNVkit](https://github.com/etal/cnvkit), [Control-FREEC](https://github.com/BoevaLab/FREEC),
[BIC-seq2](http://compbio.med.harvard.edu/BIC-seq/) and [cn.MOPS](https://bioconductor.org/packages/release/bioc/html/cn.mops.html),
WISECONDOR appeared to normalize copy number data in the most consistent way &mdash; by far. Nevertheless,
as is the case with every tool, WISECONDOR has limitations of its own: the Stouffer's z-score approach is error-prone when
dealing with large amounts of aberrations, the algorithm is extremely slow (24h) when using small bin sizes (15 kb) and
sex chromosomes are not included in the analysis. Here, I present wisecondorX, an evolved WISECONDOR that aims at dealing with
previous difficulties. Main adaptations include the additional (and consistent) analysis of the X and Y chromosomes,
a CBS-based segmentation technique and a custom plotter, resulting in overall superior results and significantly lower computing times,
allowing daily diagnostic use. WisecondorX should be applicable not only to NIPT, but also to PGD, FFPE, LQB, ... etc.

# Manual

There are three main stages for using wisecondorX:
- Converting .bam files to .npz files (both reference and test samples)
- Creating a reference (using reference .npz files)  
    - **Important notes:**
        - Reference samples should be divided in two distinct groups, one for males and one for females. This is required to correctly
        normalize the X and/or Y chromosome.  
        - When the female reference is given to the [`predict`](#stage-3-predict-cnas) function, chromosome X will be analysed;
        when on the other hand the male reference is used, chromosomes X & Y are analysed. This regardless of the gender of the test case,
        although I would **never** advice to use a female reference and a male test case, or vice versa &mdash; this because numerous Y reads
        wrongly map the X chromosome. Using a matching reference, the latter is accounted for.
        - For NIPT, exclusively a female reference should be created. This implies that for NIPT, wisecondorX is not able
        to analyse the Y chromosome. Furthermore, obtaining consistent shifts in the X chromosome is only possible when the reference
        is created using pregnancies of female fetuses only.  
        - It is of paramount importance that the reference set consists of exclusively healthy samples that originate from the same 
        sequencer, mapper, reference genome, type of material, ... etc, as the test samples. As a rule of thumb, think of
        all laboratory and in silico pre-processing steps: the more sources of bias that can be omitted, the better.  
        - Try to include at least 50 samples per reference. The more the better, yet, from 200 on it is
        unlikely to observe additional improvement concerning normalization.  
- Predicting CNAs (using the reference and test .npz cases of interest)

### Stage (1) Convert .bam to .npz

`python2 wisecondorX.py convert input.bam output.npz [-optional arguments]`  
  
<br>Optional argument<br><br> | Function
:--- | :---  
`-binsize x` | Size per bin in bp, the reference bin size should be a multiple of this value (default: x=5e3)  
`-retdist x` | Max amount of bp's between reads to consider them part of the same tower (default: x=4)  
`-retthres x` | Threshold for a group of reads to be considered a tower. These will be removed (default: x=4)  

&rarr; Bash recipe (example for NIPT) at `./pipeline/convert.sh`

##### Alternatively, convert (old) WISECONDOR .npz to wisecondorX .npz

`python2 wisecondorX.py reformat input.npz output.npz`

### Stage (2) Create reference

`python2 wisecondorX.py newref reference_input_dir/*.npz reference_output.npz [-optional arguments]`  
  
<br>Optional argument<br><br> | Function
:--- | :---  
`-gender x` | The gender of the samples at the `reference_input_dir`, female (F) or male (M) (default: x=F)  
`-binsize x` | Size per bin in bp, defines the resolution of the output (default: x=1e5)  
`-refsize x` | Amount of reference locations per target (default: x=300)  
`-cpus x` | Number of threads requested (default: x=1)  

&rarr; Bash recipe (example for NIPT) at `./pipeline/newref.sh`

##### When the gender is not known, wisecondorX can predict it

`python2 wisecondorX.py gender input.npz [-optional arguments]`  

<br>Optional argument &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Function
:--- | :---  
`-cutoff x` | Y-read permille cut-off: above is male, below is female. Note that for NIPT, this does not allow to distinguish male from female fetuses (default: x=2.5; optimized for Bowtie2 hg38 mapping)  

### Stage (3) Predict CNAs  

`python2 wisecondorX.py predict test_input.npz reference_input.npz output_id [-optional arguments]`  
  
<br>Optional argument &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Function  
:--- | :---  
`-minrefbins x` | Minimum amount of sensible reference bins per target bin (default: x=150)  
`-maskrepeats x` | Regions with distances > mean + sd * 3 in the reference will be masked, number of masking cycles (default: x=5)  
`-alpha x` | P-value cut-off for calling a CBS breakpoint (default: x=1e-4)  
`-beta x` | Number between 0 and 0.5, defines the trade-off between sensitivity and specificity for aberration calling. If e.g. beta=0.1, ratios between 0.95 and 1.05 (1.05 - 0.95 = 0.1) are seen as non-aberrant (default: x=0.05)  
`-blacklist x` | Blacklist that masks additional regions in output, requires header-less .bed file. This is particularly useful when the reference set is a too small to recognize some obvious regions (such as centromeres; example at `./blacklist/centromere.hg38.txt`) (default: x=None)  
`-bed` | Outputs tab-delimited .bed files (trisomy 21 NIPT example at `./example.bed`), containing all necessary information  **(\*)**
`-plot` | Outputs custom .png plots (trisomy 21 NIPT example at `./example.plot`), directly interpretable  **(\*)**  

<sup>**(\*)** At least one of these output formats should be selected</sup>  

&rarr; Bash recipe (example for NIPT) at `./pipeline/predict.sh`

# Parameters

The default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x depth; sWGS) and reference bin sizes 
ranging from 50 to 200 kb. When increasing the reference bin size (`-binsize`), I recommend lowering the reference locations 
per target (`-refsize`) and the minimum amount of sensible reference bins per target bin (`-minrefbins`). Further note that a
reference bin size lower than 15 kb is not advisable, unless a higher sequencing depth was used.  
**Important note**  
Concerning the vast majority of applications, the `-alpha` parameter should not be tweaked. The `-beta` parameter on the contrary
should depend on your type of analysis. For NIPT, its default value should be fine. However, for gDNA, when mosaicisms are of no interest,
it could be increased to its maximum, being 0.5. When the fetal (NIPT) or tumor (fresh material, FFPE, ...) fraction is known, this parameter is optimally
close to the fraction / 2. If you have any doubts about this argument, a default `-beta` should still be fine when a good and large reference set was created,
irrespective the type of analysis.  

# Underlying algorithm

To understand the underlying algorithm, I highly recommend reading [Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809); and its
update shortly introduced in [Huijsdens-van Amsterdam et al (2018)](https://www.nature.com/articles/gim201832.epdf).
Some adaptations to this algorithm have been made, e.g. additional median centering and variance stabilization (log2) on final ratios, removal of
less useful plot and Stouffer's z-score codes, addition of the X and Y chromosomes, inclusion of CBS and plot codes, and &mdash; last but not least &mdash;
restrictions on within-sample referencing:  

![Alt text](./figures/within-sample-normalization.png?raw=true "Within-sample normalization in wisecondorX")

# Dependencies

- R (v3.4) packages
    - jsonlite (v1.5)
    - png (v0.1-7)
- R Bioconductor (v3.5) packages
    - DNAcopy (v1.50.1)
- Python (v2.7) libraries
	- scipy (v1.0.0)
    - scikit-learn (v0.19.1)
    - pysam (v0.13)
    - numpy (v1.13.3)

And of course, other versions are very likely to work as well.  