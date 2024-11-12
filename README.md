# perffseq_reproducibility
Analyses to reproduce PERFF-seq manuscript (Abay _et al. Nature Genetics_ 2024).

## Upstream processing
Fastq data were processed using default parameters with [CellRanger v7.2](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in). 
Using a new-ish version (v7+) is necessary to get the half-mapped summary statistic, which was 
used extensively for the technology development inferences. 

Raw .fastq data is [available on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262355)
that can serve as input into CellRanger for reproducing the processed data in this work. 

For all analyses, we used the Human or Mouse v1.0.1 probeset. All experiments
were performed on human cells except the `mouse_brain_nuclei` experiment. 

## Usage
This repository contains all custom downstream analyses for reproducing
the figures and summary statistics per figure in the PERFF-seq paper. 

In brief, the processed counts data (generally `.h5` format) are exported from the 
CellRanger processing and organized according to the downstream analysis task. 

The key contents are summarized in these folders:

```
Figure 1 - benchmark/quality_control_techdev
Figure 2 - benchmark/fourplex
Figure 3 - logic_gating
Figure 4 - tfs_pbmcs
Figure 5 - loy
Figure 6 - human_ffpe_nuclei,mouse_brain_nuclei
```

Custom code is written in R and makes extensive use of the [Seurat](https://satijalab.org/seurat/) ecosystem.

## Resources
Here is a shared repo containing code for [custom 10x probe designs](https://github.com/clareaulab/flex-custom-probes).
Though not used in this work, these design tools are useful for creating custom Flex probes
for diverse applications (e.g., GFP expression, lncRNA detection).

### Contact: 
- [Tsion](mailto:tabay@g.harvard.edu)
- [Bob](mailto:rstickel@stanford.edu)
- [Ronan](mailto:chalignr@mskcc.org)
- [Ansu](mailto:satpathy@stanford.edu)
- [Caleb](mailto:lareauc@mskcc.org)

