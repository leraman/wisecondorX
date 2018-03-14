# Introduction  
After extensively comparing different WGS-based CNA tools, such as [wisecondor](https://github.com/VUmcCGP/wisecondor),
[QDNAseq](https://github.com/ccagc/QDNAseq), [CNVkit](https://github.com/etal/cnvkit), [Control-FREEC](https://github.com/BoevaLab/FREEC),
[CNVnator](https://github.com/abyzovlab/CNVnator), [ichorCNA](https://github.com/broadinstitute/ichorCNA), ... etc,
wisecondor appeared to normalize copy number data in the most consistent way &mdash; by far. Nevertheless,
as is the case with every tool, wisecondor has limitations of its own: the Stouffer's z-score approach is error-prone when
dealing with large amounts of aberrations, the algorithm is extremely slow (24h) when using small bin sizes (15 kb) and
sex chromosomes are not included in the analysis. Here, I present wisecondorX, an evolved wisecondor that aims at dealing with
previous difficulties. Main adaptations include the additional (and consistent) analysis of the X and Y chromosomes,
a CBS-based segmentation technique and a custom plotter, resulting in overall better results and significantly lower computing times.
WisecondorX should be applicable not only to NIPT, but also PGD, FFPE, LQB, ... etc.

# Manual

There are three main stages for using wisecondorX:
- Converting .bam files to .npz files (both reference and test samples)
- Creating a reference (using reference .npz files)  
    - **Important notes:**
        - Reference samples should be divided in two distinct groups, one for males and one for females. This is required to correctly
        normalize the X and/or Y chromosome.  
        - When the female reference is given to the [`predict`](#stage-3-predict-cnas) function, chromosome X will be analysed;
        when on the other hand the male reference is used, chromosomes X & Y are analysed. This regardless of the gender of the test case,
        although I would **not** advice to use a male reference and a female test case, or vice versa &mdash; this because numerous Y reads
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

##### Alternatively, convert (old) wisecondor .npz to wisecondorX .npz

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
`-blacklist x` | Blacklist that masks additional regions in output, requires header-less .bed file. This is particularly useful when the reference set is a too small to recognize some obvious regions (such as centromeres; example at `./blacklist/centromere.hg38.txt`) (default: x=None)  
`-json` | Outputs .json file, containing all output information  **(\*)**
`-bed` | Outputs tab-delimited .bed files, containing most import information  **(\*)**
`-plot` | Outputs custom .png plots (healthy male example at `./example.plots`), directly interpretable  **(\*)**  

<sup>**(\*)** At least one of these output formats should be selected</sup>  

&rarr; Bash recipe (example for NIPT) at `./pipeline/predict.sh`

# Parameters

The default parameters are optimized for shallow whole-genome sequencing data (0.1x - 1x depth; sWGS) and reference bin sizes 
ranging from 50 to 200 kb. When increasing the reference bin size (`-binsize`), I recommend lowering the reference locations 
per target (`-refsize`) and the minimum amount of sensible reference bins per target bin (`-minrefbins`). Further note that a
reference bin size lower than 15 kb is not advisable, unless a higher sequencing depth was used.

# Calling an aberration

WisecondorX will call segments, but will not call 'true aberrations'. I believe there are too many parameters, such as tumor content 
(LQB, FFPE, ...), fetal content (NIPT), sequencing depth, number of tests, p-value cut-off, ... etc, to reliably and automatically call an aberration
(with sWGS input, at least). Users can always exploit the "cbs_calls" **(\*)** and/or the "results_r" **(\*\*)** .json keys of `wisecondorX.py predict`'s
output to program a personal solution.  
  
<sup>**(\*)** Total collection of segments structured as \[chr, start index, end index, BM-score (median of bin wise z-scores of the segment), log2-ratio\]</sup>  
<sup>**(\*\*)** Total collection of bin wise log2-ratio's</sup>  

# Underlying algorithm

To understand the underlying algorithm, I highly recommend reading [Straver et al (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24170809).
Some minor adaptations to this algorithm have been made, e.g. additional median centering and variance stabilization (log2) on final ratios, removal of
less useful plot and Stouffer's z-score codes, addition of the X and Y chromosomes, and &mdash; last but not least &mdash;
restrictions on within-sample referencing:  

![Alt text](./figures/within-sample-normalization.png?raw=true "Within-sample normalization in wisecondorX")

# Additional dependencies

- R version 3.3 or higher
- R packages
    - jsonlite (v1.5)
    - png (v0.1-7)
- R Bioconductor packages
    - DNAcopy (v1.50.1)

And of course, other versions might work as well.  