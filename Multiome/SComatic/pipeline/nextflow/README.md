# SComatic Nextflow pipeline

SComatic is a tool for identifying somatic mutations across populations of single cell data. The Nextflow pipeline `scomatic.wf` handles everything from getting the data from iRODS to providing a few different versions of filtered mutations at the end of the SComatic workflow.

## Setup

- Have job submission prowess on the farm by completing the farm course. If you haven't done this, write servicedesk to enrol.
- Have a conda variant installed on the farm. Within said conda variant, have an environment with the package `docopt` installed - this allows you to use `jsub` rather than `bsub`, offering simpler syntax.
- Make sure that your `~/.bashrc` has the line `export LSB_DEFAULT_USERGROUP="teichlab"` in it. If you just added it, disconnect from the farm and connect again.

## Input

Running the pipeline requires you to prepare two input files.

### Mappings

The mappings file is a CSV with three columns. The first is the sample ID, the second is the location of the BAM on iRODS, and the third is the donor ID. All samples for the same donor will be merged for SComatic processing. A header line is expected, but you're not forced to name the columns in the same exact manner as in the examples.

When given an iRODS folder, `ils` it and check for the BAM's name. As a general rule of thumb, the BAMs are named `possorted_genome_bam.bam`, but in the case of multiome samples the two modalities' BAMs are `atac_possorted_bam.bam` and `gex_possorted_bam.bam`. As such, this would be an example multiome ATAC mapping file:

```
ID,Mapping_iRODS,donor
BHF_F_Hea11064670_BHF_F_Hea11031823,/seq/illumina/cellranger-arc/cellranger-arc200_count_2b5837c5328f291d93f52678086d5447/atac_possorted_bam.bam,C82
BHF_F_Hea11064671_BHF_F_Hea11031824,/seq/illumina/cellranger-arc/cellranger-arc200_count_bf8b5a2b9f635a1360a0e4fdfc10c8d9/atac_possorted_bam.bam,C85
```

And this would be its GEX counterpart:

```
ID,Mapping_iRODS,donor
BHF_F_Hea11064670_BHF_F_Hea11031823,/seq/illumina/cellranger-arc/cellranger-arc200_count_2b5837c5328f291d93f52678086d5447/gex_possorted_bam.bam,C82
BHF_F_Hea11064671_BHF_F_Hea11031824,/seq/illumina/cellranger-arc/cellranger-arc200_count_bf8b5a2b9f635a1360a0e4fdfc10c8d9/gex_possorted_bam.bam,C85
```

Please make sure you have no dashes in your sample IDs in the first column. Consider replacing them with underscores if necessary. Note that the mapping file will be subset to just samples that appear in the cell type file.

In the rare event of having the mappings on the farm rather than iRODS (the best guess as to why this would be is some form of external data), there is also a non-iRODS mode of mapping loading. Prepare the mappings file like you would normally, except with a farm path to the BAM. When running the workflow, add `--location local` to the Nextflow call.

### Cell types

The cell types file is a TSV file with two columns. The first is the cell barcode, specifically formatted `SAMPLE_BARCODE-1` (the fact it's an underscore and then a dash is important), with the `SAMPLE` matching one of the IDs provided in the mappings CSV file. Please make sure all samples from the mapping file are represented. The second column is a cell type assignment. Please keep the cell types to letters, numbers and underscores if possible - any characters not fitting this definition will get converted to underscores as part of the processing, as cell types are used to spawn file names.

The header has to be `Index` and `Cell_type` or SComatic won't understand the file. An excerpt from an example file:

```
Index   Cell_type
BHF_F_Hea11064670_BHF_F_Hea11031823_TTGCTTAGTGAGACTC-1  VentricularCardiomyocytesLeftTrabeculated-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_ATAGATGCATTGTCCT-1  MesothelialEpicardialCells-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_GTACACCCATCCCTCA-1  AtrioventricularNodeCardiomyocytes-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_GGTCGGTTCTTAGGAC-1  AtrialCardiomyocytesLeft-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_GGAACAATCAAGCTTA-1  ChromaffinCells-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_CGATCCTTCTTGTCCA-1  AtrialCardiomyocytesLeft-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_GCACATTAGGGACGCA-1  EndocardialCells-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_GCTCAACCAGGCTGTT-1  AtrialCardiomyocytesCycling-C82
BHF_F_Hea11064670_BHF_F_Hea11031823_AGTAATCGTCCTAGTT-1  GreatVesselSmoothMuscleCells-C82
```

It is recommended to include the donor ID in the cell type.

### Multiplexed data

In the event of working with pooled donors, do the following:
- Identify all samples covering a given donor pool. These will be a number of overlapping multiplexed samples (e.g. donor A + donor B, and then donor A + donor C), plus any single-donor samples that contain the donors from the multiplexing. This will allow all cells for a given donor to be analysed together.
- Mark all of the samples with a new donor name (e.g. `MULTI1`, `MULTI2`, `MULTI3` etc.) in the mappings CSV.
- Be sure to include the original donor ID for each cell in the cell types TSV file. SComatic splits the cells into cell types, so `Tcell_donorA` and `Tcell_donorB` will get independent mutation profiles.

Excerpts of files with this in action - `MULTI2` is a mixture of `Donor2392` and `Donor2394`:

```bash
$ grep "MULTI2" mappings.csv 
HCA_BN_F12922484_and_HCA_BN_F12918722,/seq/illumina/cellranger-arc/cellranger-arc201_count_203069c7060f9e1e493f20a0b051ca31/gex_possorted_bam.bam,MULTI2
HCA_BN_F12922485_and_HCA_BN_F12918723,/seq/illumina/cellranger-arc/cellranger-arc201_count_7cf9099c70409c0824c6c263080325fc/gex_possorted_bam.bam,MULTI2
HCA_BN_F12922486_and_HCA_BN_F12918724,/seq/illumina/cellranger-arc/cellranger-arc201_count_ad094b8d403f5ebc05655f1099da327e/gex_possorted_bam.bam,MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725,/seq/illumina/cellranger-arc/cellranger-arc201_count_2ec8d25de2ee618f8df84ca88fe867a5/gex_possorted_bam.bam,MULTI2
```

```bash
$ grep "HCA_BN_F12922487_and_HCA_BN_F12918725" celltypes.tsv | head
HCA_BN_F12922487_and_HCA_BN_F12918725_AAACAGCCATTGCGGT-1        Capillary_endothelial_cell-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AACGCTAGTTATAGCG-1        Arterial_endothelial_cell-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AAGAATCAGGTGAAGC-1        Capillary_endothelial_cell-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AAGCCTGTCCTACCTA-1        Paraximal_mesoderm_with_endothelial_protential-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AAGCTCCCAACACTTG-1        Capillary_endothelial_cell-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AAGGAAGCAAAGCCTC-1        Venous_endothelial_cell__BBB_-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AAGTGAAGTGTTGTAG-1        Capillary_endothelial_cell-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_AATTGTGTCATGACCG-1        Capillary_endothelial_cell-Donor2392-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_ACAACAACATGCTTAG-1        Arterial_endothelial_progenitor-Donor2394-MULTI2
HCA_BN_F12922487_and_HCA_BN_F12918725_ACCATAATCAGGAAGC-1        Capillary_endothelial_cell-Donor2394-MULTI2
```

## Running the pipeline

Once you're set up, running the pipeline is very easy. Make a folder on Lustre, then create `run_scomatic.sh`:

```bash
#!/bin/bash
set -eo pipefail

#conda activate an environment which has docopt installed so jsub works
#/software/team205/bin/jsub lsf -q week -n scomaticGEX -c 1 -m 2g -l logs "bash run_scomatic_GEX.sh" | bsub

#input files, two samples from jc48, representative of formatting
MAPPINGS=/nfs/team205/kp9/nextflow/scomatic/demo/mappings-GEX.csv
CELLTYPES=/nfs/team205/kp9/nextflow/scomatic/demo/celltypes-GEX.tsv

#GEX.json and ATAC.json - parameters configured for 10X's GRCh38-2020-A reference
PARAMETERS=/nfs/team205/kp9/nextflow/scomatic/GEX.json

mkdir -p logs

/software/team205/nextflow run /nfs/team205/kp9/nextflow/scomatic/scomatic.nf \
    --mappings $MAPPINGS \
    --celltypes $CELLTYPES \
    --projectDir $PWD \
    -params-file $PARAMETERS \
    -c /nfs/team205/kp9/nextflow/scomatic/LSF.config \
    -resume
```

Point `MAPPINGS` and `CELLTYPES` to the files you prepared and you should be able to submit the job. I tend to use `jsub` because of its simpler syntax, and `/software/team205/bin/jsub lsf -q week -n scomaticGEX -c 1 -m 2g -l logs "bash run_scomatic_GEX.sh" | bsub` should work as stated in a comment in the script.

### GEX/ATAC

The pipeline needs slightly different input for GEX/ATAC, and there are two parameter files in `/nfs/team205/kp9/nextflow/scomatic/` calibrated for the two modalities. The ATAC mode is calibrated for multiome ATAC. If processing multiome samples, run the pipeline in GEX and ATAC modes separately, then join the two sets of results in post-processing.

### Older references

The pipeline's parameter files are calibrated for 10X's GRCh38 2020-A reference. If you had your data mapped to an older reference (e.g. 10X's GRCh38 3.0.0), do the following:
- Copy the relevant modality JSON file, as you'll need to edit it.
- Edit the `"editing"` (GEX only) and `"pons"` paths, replacing the `.txt` at the end with `.pre2020.txt`. Similarly, replace the `.bed` with `.pre2020.bed` in the `"bed"` path.
- Replace `"genome"` with a path to the genome fasta file from your reference. Please make sure it's accompanied by a `.fai` file in the same folder (which can be created via `samtools index` if necessary).
Change the `PARAMETERS` path in the script to your edited JSON file.

## Output

The pipeline yields the following output formatting:

```
C82-GEX/
├── C82.calling.step2.intersect.tsv
├── C82.calling.step2.pass.tsv
├── C82.calling.step2.tsv
├── C82.coverage_cell_count.per_chromosome.report.tsv
├── C82.coverage_cell_count.report.tsv
└── cell_callable_sites
    ├── C82.AtrialCardiomyocytesCycling_C82.SitesPerCell.tsv
    ├── C82.AtrialCardiomyocytesLeft_C82.SitesPerCell.tsv
    ├── C82.AtrialCardiomyocytesRight_C82.SitesPerCell.tsv
    ├── C82.AtrioventricularNodeCardiomyocytes_C82.SitesPerCell.tsv
    ├── [...]
```

`calling.step2.tsv` is the output of the last step of SComatic. In [the example](https://github.com/cortes-ciriano-lab/SComatic/blob/main/docs/SComaticExample.md), the author recommends subsequently intersecting this file with a BED of high quality regions of the human genome, and then filtering to just the mutations that `PASS` all of the filters. These two sets of mutations are provided in `calling.step2.intersect.tsv` and `calling.step2.pass.tsv` respectively. James consulted the exact use scenario of the output with the author, and starts from `calling.step2.intersect.tsv`, doing a less restrictive filtering.

Once your SComatic run finishes, please delete the `work` directory that gets created. While Nextflow kind of cleans up after itself, its power to do so is limited to the exact processes that were ran as part of that execution. For example, if you needed to restart, the cached jobs from earlier will have their various files left untouched. Also even when cleaning up, Nextflow leaves behind the various command files, creating a significant number of files that mean very little to you at that stage.

## Single cell genotyping

Once you have the output of SComatic, you can optionally identify genotypes at the single cell level. This is not done by default as the output is far bigger than the normal SComatic output (usually at least 10x as large, possibly even more!). As such, **it is recommended to subset your cell types file to just cell types of interest before proceeding.**

The requisite bash script is very similar to the normal SComatic one, but requires you to include a path to the folder where you originally ran SComatic, the folder with the various `DONOR-MODALITY` subfolders that got created, as `MUTATIONS`. The pipeline will then load the files it needs from there.

An example `run_scomatic_GEX_genotypes.sh`:

```bash
#!/bin/bash
set -eo pipefail

#conda activate an environment which has docopt installed so jsub works
#/software/team205/bin/jsub lsf -q week -n scomaticGEXgenotypes -c 1 -m 2g -l logs "bash run_scomatic_GEX_genotypes.sh" | bsub

#input files, two samples from jc48, representative of formatting
MAPPINGS=/nfs/team205/kp9/nextflow/scomatic/demo/mappings-GEX.csv
CELLTYPES=/nfs/team205/kp9/nextflow/scomatic/demo/celltypes-GEX.tsv
#scomatic pipeline output for mutations to use
MUTATIONS=/nfs/team205/kp9/nextflow/scomatic/demo

#GEX.json and ATAC.json - parameters configured for 10X's GRCh38-2020-A reference
PARAMETERS=/nfs/team205/kp9/nextflow/scomatic/GEX.json

mkdir -p logs

/software/team205/nextflow run /nfs/team205/kp9/nextflow/scomatic/scomatic.nf \
    -entry genotypes \
    --mappings $MAPPINGS \
    --celltypes $CELLTYPES \
    --mutations $MUTATIONS \
    --projectDir $PWD \
    -params-file $PARAMETERS \
    -c /nfs/team205/kp9/nextflow/scomatic/LSF.config \
    -resume
```

## Scripts

Just in case, the relevant scripts and configuration files are included in this folder. You shouldn't need them though, as they live in `/nfs/team205/kp9/nextflow/scomatic/` on the farm and should be accessible. There are a couple changes relative to stock SComatic, captured in [this pull request](https://github.com/cortes-ciriano-lab/SComatic/pull/9).