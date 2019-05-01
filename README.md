# usl-count #

Gene location-aware accurate quantitation of gene expression from
high-throuput sequencing alignments.

# Installation #

    pip install -e git+https://github.com/mikpom/uslcount.git#egg=uslcount
    
`--prefix` option of pip install can be used to specify custom
installation root directory.

Add `--upgrade` to upgrade existing installation. 
    
# Intended use #

Package can be used to quantitate reads by gene features from BAM
alignments.  For unstranded libraries on top of gene expression
characteristics package can give the confidence score idicating how
likely the counts are compromised (overestimated) due to expression
from the opposite strand. One might need to manually inspect those
with low confidence score and maybe consider using splice junction
counts (also contained in output file) as an expression estimate of
these genes.

# Typical workflow #

It is assumed that NGS reads were already aligned to the genome using
aligner of choice (e.g. STAR).  Prior to running quantitation analysis
a genomic database should be created which will be then used for
processing of BAM files.

Available tasks can be listed by typing:

    python3 -m uslcount
    
## Building genomic DB ##
    
To build genomic database run 

    python3 -m uslcount build /path/to/GTF /path/to/GenomicDataDir
    
## Counting reads from BAM ##

To count reads run

    python3 -m uslcount count -out outfile GenomicDataDir BAM

where GenomicDataDir indicated path to data directory created at the `build` step. 

To see other options

    python3 -m uslcount count -h
    
#### Effective gene length ####

To convert count values to TPM one needs to know gene length parameter
to make comparisons between different genes more adequate. Actual span
from which counts are gathered is dependent on whether the analysis
was stranded or not. If counting was unstranded then corresponding
gene lengths can be obtained from the file `len_ust_unq`. For stranded
analysis lengths are stored in the file `len_st_unq`. The latter
lengths are equal or larger than the former. Both files can be found
in genomic directory indicated at the `build` step.

Note that for some genes length equals to zero meaning there is no
non-overlapping exons to distinguish that gene from the other. In this
case count is always zero.

## Analyzing BAM file obtained from reads of unstranded library ##

To get confidence score of counts obtained from unstranded libraries
`analyze` task should be invoked. Usage is very similar to count task:

    python3 -m uslcount analyze -out outfile GenomicDataDir BAM

This will output counts, junction counts and a confidence score (int
from 1 to 10) indicating how robust expression value is. 10 means
expression value is reliable, 1 means -- least reliable. Junction
counts can be used as a proxy for gene count in differential
expression analysis in those cases when gene's expression biased by
opposite strand transcription.
