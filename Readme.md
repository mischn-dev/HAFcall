# HAFcall.jl - Haplotype Analysis of low coverage sequencing data

The package performs a frequency calling strategy for pooled population / offspring samples for low coverage sequencing data

Analyzing Population structures and estimate the allele frequency can be extremely expensive or very limited in its statement.
This package highlights an approach to estimate the haplotype frequency with very high security on an ultra-low cost approach using pooled sampling of hundreds of genotypes together. 

## Resources 
 
- **Publication** <http://linktopub>

Three things are necessary to implement
- a **vcf** file including the polymorphisms for all samples tested
- a **overview** file where all information about the samples is stored 
- a **reference** file derived from a database, like a GFF3 file for genes 

## Installation

Within julia, execute
```julia
using Pkg; Pkg.add("HAFcall"); using HAFcall
```

### Troubleshooting

If an error occurs, try `Pkg.resolve()` and repeat `using HAFcall`


## Purpose of this Package

The package has a couple of functions that need to be chained.

The **first step** `AFcall` extracts the SNP information for all the samples from the Vcf file and stores these in a Dictionary - one entry for each sample
The output is stored as DataFrames and can be accessed by `DictName[:Set1]`, where `:Set1` is equivalent to the sample name if the **overview** file.

The **second step** is to `locate` the Polymorphisms. The function compares the position of a given object, like a gene (derived from e.g. a **GFF3** file), and the polymorphisms and create a dictionary entry label `:Location`, where each polymorphism is tested for its link to an object. 

The **third step** performs the `haplotyping`, where the information of multiple polymorphisms is melted together to one haplotype frequency value for each object. 

Merging the allele frequency of multiple low coverage SNPs / INDELS located in the sample genomic object that is expected to be linked (linkage disequilibrium), will result in a much higher coverage and therefore much better estimate of the frequency [a gene should have a LD value of 1, crossing overs are unlikely because the region is very small]
To increase the power of this method, the genomic object can be extended into the nearby genomic region. Due to linkage, these Polymorphism should still represent the same block and can be utilized to increase the amount of Polymorphism and coverage to receive an even better haplotype frequency estimate.


# Usage

## required packages:

these should be already loaded, but if there are any problems, try to install the following packages first:
```julia
] add  DelimitedFiles DataFrames Statistics CSV ProgressMeter Crayons StatsBase
```

The package supports multi-threading, so make sure you have loaded julia in multicore mode - check with `Threads.nthreads()` , if there is only 1 core active run 

	JULIA_NUM_THREADS=8 julia or julia -t 8

where 8 is a dummy for the number of cores. Set it to a value your machine is running on (e.g., an early Intel i7 has 8 threads, an i5 only 4)
Alternative, if you run julia in **Atom**, check out the **julia-client** options
Also, make sure **julia** can be found by your system. If you saved **julia** in the Documents folder, add the path to it like `C:/Users/me/Documents/julia/bin/julia.exe`

## required Data

### Overview file

the functions need to know which column in the vcf dataset belongs to which data. Furthermore, the environment and the replicates have to be specified for some functions.

Example: 
| Name | Env  | Generation  |
| ------- | --- | --- |
| Scarlett | 0 | 0 |
| Trabant | 0 | 0 |
| F3R1 | 1 | 3 |
| F3R2 | 1 | 3 |
| F3R3 | 1 | 3 |


The Parents / Founders / near Relatives should be coded both 0 for Env and Generation. Env stands for environment. 
When Env and Generation information is derived, additional statistical analysis can be performed.

The file has to be *TAB separated* e.g. a text file (.txt)

### Reference file

A reference file is needed to assign the Polymorphism to a genomic Object, like a Gene, Marker, Exon, UTR, etc. the file should have at least the following information:
- Chromosome
- Start
- End

Example 1 *(presplit = false in function location)*:

| Gene_ID   | chromosome:start-stop  |
| ------------- | ------------- |
| HORVU1Hr1G000010  | chr1H:41961-45310   |
| HORVU1Hr1G000020  | chr1H:47115-50750    |
| HORVU1Hr1G000030  | chr1H:50869-52927    |

Example 2 *(presplit = true in function location)* : 

| Marker | Chr  | start  | end |
| :-----: | :-: | :-: | :-: |
| BOPA2_12_10420 | chr1H  | 20000  | 71484 |
| BOPA1_3107-422 | chr1H  | 71563   | 151117   |
| SCRI_RS_204276 | chr1H  | 256663  | 267591 |
| BOPA2_12_30653 | chr1H  | 273622  | 366508 |

additional information can be added in more columns, like:

| Marker | Chr  | start  | end | Pos_genetic | description | GO-IDs
| :-----: | :-: | :-: | :-: | :-: | :-: | :-: |
| BOPA2_12_10420 | chr1H  | 20000  | 71484 | 0.102 | ribosomal protein | GO:000347
| BOPA1_3107-422 | chr1H  | 71563   | 151117   | 0.132 | HMA | GO:004547
| SCRI_RS_204276 | chr1H  | 256663  | 267591 | 0.212 | Calmodulin | GO:0007874
| BOPA2_12_30653 | chr1H  | 273622  | 366508 | 0.845 | NAC | GO:0004748


The start and end position do not need to differ, but the file has to be *TAB separated \t*

### VCF file

the vcffile can be produced by various SNP calling software tools. This tool was designed based on `SAMtools` and `bcftools` Polymorphism calling using *mpileup* and *call*
since version 1.9, both *mpileup* and *call* have shifted to bcftools

To receive the correct data format of the vcf file, us code that results in the same output as *Samtools / Bcftools version 1.8* [-t AD option is crucial]
```shell
samtools mpileup -uvf Ref_dna.toplevel.fa -q 25 -Q 30 -t AD first.bam second.bam third.bam fourth.bam | 
bcftools call -vm -Ov > Polymorphismfile.vcf 
```

The calling output should look like this for each Polymorphism detected:

| Chr | Pos  | #3  | Ref | Alt | Qual | #7 | Info | Style | Sample 1 | Sample 2  | Sample 3 | Sample 4
| :-----: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-: | :-:| :-:|
|chr1H|	254597	| . | A	| G	| 999 | . |	DP=220;VDB=0.001[..];DP4=75,56,23,30;MQ=59 |	GT:PL:AD |	0/0:0,33,255:11,0	|1/1:152,15,0:0,5	| 0/1:76,0,171:7,3 | 0/1:10,0,203:8,1

the file does not have to be converted or transformed in any way when the AD format is matched
The table above only is used to highlight the required format - the actual file should be a *raw vcf file* as generated by the vcf caller.  

## Workflow

The package has several functions, some of them have to be called in a row to receive a final result. 

### First function

 always is `Allelefrequcalling` - this performs the conversion of the vcf file into Dataframes, separate for each tested sample. 

```julia
Allelfreqcalling(vcffile, genolist, minreaddepth=1, minqual=30, maxqual=1000,posE1=1, posE2=2) \n

overviewfile = "info.txt"
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2)
```

some adjustments can be made, like minimum coverage or quality. Additionally, it has to be specified in which row in the overviewfile the parents can be found

the function will output a *Dictionary* with entries - each sample has an own data frame in the Dictionary that can be accessed by 
```julia
freq[:F3R1]
```

### Second function 

always is `locate` - the genomic objects will be linked to the Polymorphisms 

```julia
locate(info, freqfile, genolist, presplit=false,  extention=true, gap = 0.45) 

freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2)
locate(reference_file, freq, overviewfile)
```

- presplit - false => reference file looks like *Example 1*
           - true  => reference file is already split for the position, like in *Example 2*
- extension => shall the Objects be wide up or shall they keep their start and end position - true if extension should be performed (recommended), false if not 
  -> including the extension will give better results for the haplotype estimate 
- gap => how far should the Object be extended? e.g., how far a gene should be extended into the intergenic region (default = 0.45 - covering 90% of the intergenic region for up- and downstream gene) 
 -> do not use more than 0.5 -> the regions will overlap (not recommended, but can be done) 

a dictionary entry *Location* is added to the previously created Dict() by *AFcall*


### following functions

the `HAFcall` package has some more functions:

- haplotyping
- contig
- readDict
- writeDict
- stacker
- annotate

#### Haplotyping

`haplotyping` creates haplotypes of the polymorphisms using the `locate` generated genomic object information

 haplotype have been proven to give a higher consistency in the calling of the correct allele frequency 
 for the haplotype calculation 2 ways can be applied, the first is referring to the read depth given, so higher rd will contribute more to the Frequency (weighted = true - default, recommended) 
 second way weights all polym. equally, not depending if they are unweighted in read depth (weighted = false)

 - gene = haplotyping by Genes (Column in Location Dict entry has to have the name "Gene_ID") 
 - Marker = haplotyping by Marker (Column has to have the name "Marker")
 - GO = haplotyping by go terms (Go Enrichment) (column has to have the name "GO-IDs") 
 - PF = haplotyping by PF terms (column has to have the name "PFAM-IDs") 

 *haplotyping(freqfile,genolist,haplotype="gene", weighted=true)*

```julia
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2) 
locate(reference_file, freq, overviewfile) 
haplotyping(freq, overviewfile, "gene")
```

The output is stored as additional Dictionary entries, that do have the ending *_Haplotypes*

```julia
freq[:F3R1_Haplotypes]
``` 

`contig` performs haplotype creation based on a window approach - this can be used if genomic objects are not given by any reference or other source.

if no information on the gene or marker position can be retained from any source, the haplotypes can be created by generating contigs according to a certain size. 
inbreeding species can use bigger contigs, while outbreeding species should have smaller contig size 

default contig size = 10 mio bp 

```julia
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2) 
contig(freq, overviewfile, 10000)
``` 

The output is stored as additional Dictionary entries, that do have the ending *_Contigs*

```julia
freq[:F3R1_Contigs]
``` 

#### Read / write 

with `readDict` and `writeDict`, Frequency calls generated by the previous mentioned function can be stored to disk or read from disk.

```julia
overviewfile = "info.txt"
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2)

writeDict("C:/Users/me/Desktop/folder", freq)

readDict("C:/Users/me/Desktop/folder")

```

the `readDict` function can also only read files that do have a specific character of interest in their name - e.g. `readDict("C:/Users/me/Desktop/folder", "F")` to read files that do have a capital ?F? in its name


#### Helper functions

`stacker` can merge the replicates together

check the Dict. and check how the naming is done - can be either _Haplotypes (default) or _Contigs 
you can merge replicates if 2 or 3 replicate sets are present - for more replicates please check the code and adjust by yourself 

```julia
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2) 
locate(reference_file, freq, overviewfile)
haplotyping(freq, overviewfile, "gene")
stacker(freq, overviewfile, "_Haplotypes")
```


`annotate` adds the position information to the Haplotype DataFrames 


```julia
freq  = AFcall(SNPfile.vcf, overviewfile, 3,30,1000,1,2) 
locate(reference_file, freq, overviewfile) 
haplotyping(freq, overviewfile, "gene") 
annotate( freq[:F2E2_Haplotypes], freq, "Gene_ID")

```
