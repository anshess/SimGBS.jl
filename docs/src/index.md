# SimGBS

A Julia package for simulating Genotyping-by-Sequencing (GBS) data. 


## Background 


Genotyping-by-Sequencing (GBS) offers a cost-effective way to obtain genomic information of many samples, it has been applied in various research and commercial programmes. Several downstream analytical methods have been developed to handle GBS data. 

However, there is limited resource available to 1. assess different GBS experimental designs and 2. evaluate exisiting bioinformatics workflows and statistical methods designed for GBS data. This can affect the reproducibility of any proposed GBS study. 

A method to simulate GBS data is therefore essential to guide future applications of GBS. Here we present SimGBS, a simple yet versatile method to model GBS process and simulate GBS data. 



## Overview

`SimGBS` can be implemented in three simple steps 

 - *Genome*: generate GBS fragments via _in slico_ digestion using common restriction enzyme(s) 
 - *Population*: define population structure by implementing the _gene-drop_ method  
 - *Data*: produce GBS reads under different sequencing efforts, following the systematic sampling of genetic variants 


## Getting Started 

### Installation 
```{julia}
julia> Pkg.add("SimGBS")
```


### Examples

```@docs
GBS(totalQTL::Int64, totalSNP::Int64, muSNPdensity::Float64, sigmasqSNPdensity::Float64, winSize::Int64, muAlleleFreq::Float64,sigmasqAlleleFreq::Float64, re, meanDepth::Float64,   barcodeFile::String, plotOutput::Bool, writeOutput::Bool,outputOnlyGBS::Bool)
```



## Arguments
* `totalQTL`: The number of QTL to be simulated
* `totalSNP`: The number of SNPs to be simulated (set to "0" if sampling SNP positions based on density)
* `muDensity`: Average SNP density
* `sigmasqDensity`: Variance of SNP density
* `winSize`: Size of window and bin for sampling SNP positions
* `muAlleleFreq`: Average allele frequency
* `sigmasqAlleleFreq`: Variance of allele frequency
* `re`: Restriction Enzyme to be used
* `barcodeFile`: File contains GBS barcodes
* `plotOutput`: Logical, TRUE if graphical outputs are required
* `writeOutput`: Logical, TRUE if text outputs are required