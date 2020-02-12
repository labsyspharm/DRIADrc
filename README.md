# DRIADfigs: scripts to fully reproduce manuscript figures

This repository is a companion to Drug Repurposing In Alzheimer's Disease (DRIAD), a machine learning framework for quantifying potential associations between drugs and Alzheimer's Disease.

Create a [synapse.org](https://www.synapse.org/) account, if you don't have one already. Ensure that you have LaTeX and R  installed. To reproduce figures from the manuscript, run the following from the command line:

``` sh
Rscript setup.R
make
```

This will generate `figures.pdf` and `supplement.pdf`.

## Funding

We gratefully acknowledge support by NIA grant *R01 AG058063: Harnessing Diverse BioInformatic Approaches to Repurpose Drugs for Alzheimers Disease.*


