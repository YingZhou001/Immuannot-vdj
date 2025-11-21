# Immuannot-vdj pipeline introduction

Immuannot-vdj is designed to analyze Immunoglobulin (IG) and T cell receptor (TR) genes.
It can be used to annotate gene structures such as leader parts and RSS for assemblies and long reads (HiFi, Ont R10, RNA transcripts).

A manuscript is in preparing.

Last update date : 11/21/2025

Content
--------

- [Detection Strategy](#detection-strategy)
- [Installation](#installation)
    - [Download](#download)
    - [Dependency](#dependency)
- [Input and Output](#input-and-output)
- [An Example](#an-example)
- [Limitation](#limitation)
- [Todo list](#todo-list)

# Detection Strategy

The major difficulty for annotating V, D, and J coding regions from whole genome
assembly is the short sizes of those genes, especially for the D coding region which is about ~10-20bp. 
One solution is to perform basepair-to-basepair alignment to the restricted region (
extracted from whole genome seqeunce), together with careful evaluation including the flanking sequences
to get the most reliable annotation on those coding segments.[1]
The other way is to extend the coding seqeunces with flanking RSS regions so
that extra ~50-200bp will be added which directly increases the alignment
specitivity.
Immuannot-vdj employs the second way, by generating the V,D,J reference
sequences from IPD-IMGT/LIGM-DB, it can achieve high efficiency and accuracy at the
same time when annotating human assembly.

| Segment | Length | +flanking area |
| --- | --- | --- |
| V | ~300bp | ~500bp |
| D | ~10-20bp | ~80bp |
| J | ~50bp | ~100bp |

For indetified joint region, we suggest to use TRUST4's annotator to resolve the CDR3 sequences. [2]

[1] Lees, W.D., Saha, S., Yaari, G. and Watson, C.T., 2024. Digger: directed annotation of immunoglobulin and T cell receptor V, D, and J gene sequences and assemblies. Bioinformatics, 40(3), p.btae144. \
[2] Song, L., Cohen, D., Ouyang, Z., Cao, Y., Hu, X. and Liu, X.S., 2021. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nature methods, 18(6), pp.627-630.

[\[top\]](#content)

# Installation

## Download

1) Download Immuannot-bt:

```bash
git clone https://github.com/YingZhou001/Immuannot-vdj.git
```

2) Download IMGT reference data set

```bash
wget https://ftp.ebi.ac.uk/pub/databases/imgt/ligm/imgt.dat.Z -O dat.Z 
```

3) Prepare the reference data set for Immuannot-bt

```bash
scriptdir=Immuannot-bt/script
refpref=ref
dat=dat.Z
bash ${scriptdir}/mkref.sh ${dat} ${refpref}
## reference will be outputing into a folder named as ${refpref}.${date}
```

[\[top\]](#content)

## Dependency

Immuannot-bt depends on:

* minimap2 (2.27-r1193) 
* seqtk (1.4-r130-dirty)

All source codes and pre-build binaries can be found in folder 'dep-tools'. 
If any pre-build binary does work, you may need to re-build it from the source
and copy it to the folder 'dep-tools/bin'.

[\[top\]](#content)

## Options

```bash
  Usage: bash ${SCRIPTPATH}/immuannot-vdj.sh  [OPTIONS] value
      [ --asm/--hifidna/--rna/--ontuldna   flag indicating assembly or hifi reads as input     ]
      [ -r | --refdir <string> references directory                             ]
      [ -t | --target <string> target input (.fa, fa.gz, .fq, .fq.gz)           ]
      [ -o | --outdir <string> output directory (optional, default ./output)    ]
      [ -n | --nthread <int>   num of thread (optional, default 3)              ]
```

# Input and output

Immuannot-vdj takes fasta or fastq as input, and generates a gff file for
assembly-based annotation

[\[top\]](#content)

# Examples

```bash
cd test-pub
script=../src/immuannot-vdj.sh # make use this path is correct
refdir=ref.2025-11-21 # make sure the reference folder is correct

## run on hifi dna reads
  tar=dna.fa.gz
  outdir=hifi_dna-test
  nthread=10
  time bash ${script} --hifidna -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}


  ## run ont UL dna reads
  tar=ul_dna.fa.gz
  outdir=ul_dna-test
  nthread=10
  time bash ${script} --ontuldna -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}

  ## run rna transcripts reads
  tar=rna.fa.gz
  outdir=rna-test
  nthread=10
  time bash ${script} --rna -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}

  ## run on asssembly
  tar=asm2.fa.gz
  outdir=asm-test
  nthread=10
  time bash ${script} --asm -r ${refdir} -t ${tar} -o ${outdir} -n ${nthread}

```

[\[top\]](#content)

# Limitation

[\[top\]](#content)
