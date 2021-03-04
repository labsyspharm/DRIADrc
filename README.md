# DRIADrc: full reproducibility of manuscript figures

This is a companion to Drug Repurposing In Alzheimer's Disease (DRIAD), a machine learning framework for quantifying potential associations between drugs and Alzheimer's Disease.

* Nature communications article: https://www.nature.com/articles/s41467-021-21330-0
* Primary DRIAD repository: https://github.com/labsyspharm/DRIAD
* DRIAD as a webapp: https://labsyspharm.shinyapps.io/DRIAD/

## Reproducing the figures

All scripts and data needed to fully reproduce the figures are encapsulated inside a Docker container. Ensure that you have [Docker installed](https://docs.docker.com/get-docker/). Then run the following command:

```
docker run -v "$PWD":/output/ --rm labsyspharm/driadrc make
```

This will generate `figures.pdf` and `supplement.pdf` in your current directory.

## Funding

We gratefully acknowledge support by NIA grant *R01 AG058063: Harnessing Diverse BioInformatic Approaches to Repurpose Drugs for Alzheimers Disease.*


