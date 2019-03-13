# awesome-CRISPR
List of softwares/websites/databases for genome engineering, including (but not limited to) guide design, genome editing outcome, screening analysis, etc. 

This collection is inspired by [awesome-single-cell](https://github.com/seandavi/awesome-single-cell).

## Contents

- Guide design
  - Screening library design
- Genome editing outcomes and predictions
- Screening analysis

## Guide design

- [CRISPR-DO](http://cistrome.org/crispr/) - [webserver] - A web application for the design and optimization of guide sequences that target both coding and non-coding regions in spCas9 CRISPR system across human, mouse, zebrafish, fly and worm genomes.
- [SSC](http://cistrome.org/SSC/) - [webserver] - A sequence model for predicting sgRNA efficiency in CRISPR/Cas9 knockout experiments.

### Screening library design

- [CRISPR-FOCUS](http://cistrome.org/crispr-focus/) - [webserver] -  A web-based platform to search and prioritize sgRNAs for CRISPR screening experiments. 
- [GPP Web Portal](https://portals.broadinstitute.org/gpp/public/) - [webserver] -  A web-based platform to generate matching sgRNA knockout (CRISPRko) designs for a mouse or human gene, transcript or target sequence.
- [pgRNADesign](https://bitbucket.org/liulab/pgrnadesign.git) - [python] -  A algorithm to design paired gRNAs for knocking out long non-coding RNAs (lncRNAs).

## Genome editing outcomes and predictions

- [CRISPResso](https://github.com/lucapinello/CRISPResso) - [python, webserver] - A software pipeline for the analysis of targeted CRISPR-Cas9 sequencing data. This algorithm allows for the quantification of both non-homologous end joining (NHEJ) and homologous directed repair (HDR) occurrences.

## Screening analysis

- [JACKS](https://github.com/felicityallen/JACKS) [python] - A Bayesian method that jointly analyses screens performed with the same guide RNA library.
- [MAGeCK](https://bitbucket.org/liulab/mageck) - [python] - Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (MAGeCK) for prioritizing single-guide RNAs, genes and pathways in genome-scale CRISPR/Cas9 knockout screens. 
- [MAGeCK-VISPR](https://bitbucket.org/liulab/mageck-vispr) - [python] - A comprehensive quality control, analysis and visualization workflow for CRISPR/Cas9 screens.
- [MAGeCKFlute](https://bitbucket.org/liulab/mageckflute/) - [R] - A pipeline for perform computational analysis of CRISPR screens. MAGeCKFlute combines the MAGeCK and MAGeCK-VISPR algorithms and incorporates additional downstream analysis functionalities.
- [STARS](https://portals.broadinstitute.org/gpp/public/software/stars) - [python] - A gene-ranking algorithm for genetic perturbation screens, computing a score for genes using the probability mass function of a binomial distribution. To analyze either shRNA or sgRNA based screening data.
- [BAGEL](https://sourceforge.net/projects/bagel-for-knockout-screens/) - [python] - A algorithm is designed to identify negatively selected genes, by calculating a Bayes factor for each gene representing a confidence measure that the gene knockout results in a fitness defect. Bayesian analysis of gene knockout screens using pooled library CRISPR or RNAi.
- [casTLE](https://bitbucket.org/dmorgens/castle) - [python] - Based on empirical Bayesian framework to account for multiple sources of variability, including reagent efficacy and off-target effects.
- [CERES](https://depmap.org/ceres/) - [R] -  A algorithm to estimate gene-dependency levels from CRISPR-Cas9 essentiality screens while accounting for this effect.
