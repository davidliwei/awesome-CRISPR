# awesome-CRISPR
List of softwares/websites/databases for genome engineering, including (but not limited to) guide design, genome editing outcome, screening analysis, etc. 

This collection is inspired by [awesome-single-cell](https://github.com/seandavi/awesome-single-cell).

## Contents

- Guide design
- Genome editing outcomes and predictions
- Screening analysis

## Guide design

- [CasFinder](http://arep.med.harvard.edu/CasFinder/) - [python] - A algorithm for identifying specific Cas9 targets in genomes.
- [CCtop](https://crispr.cos.uni-heidelberg.de) - [webserver] - A algorithm to predict CRISPR/Cas9 target.
- [CHOPCHOP](http://chopchop.cbu.uib.no/index.php) - [webserver] - A web tool for selecting target sites for CRISPR/Cas9, CRISPR/Cpf1.
- [CRISPR Library Designer](https://github.com/boutroslab/cld_docker) - [Software] - A software for the multispecies design of sgRNA libraries.
- [CRISPR MultiTargeter](http://www.multicrispr.net/index.html)- [webserver]a web-based tool for automatic searches of CRISPR guide RNA targets. It can find highly similar or identical target sites in multiple genes or transcripts or design targets unique to particular genes or transcripts.
- [CRISPR RGEN Tools](http://www.rgenome.net) - [webserver] - A algorithm can dentify of RGEN off-target sites without limitation by the number of mismatches and allows variations in PAM sequences recognized by Cas9. Meanwhile, it can search for RGEN targets with low potential off-target effects and high knock-out efficiencies in exon.
- [CRISPRdirect](http://crispr.dbcls.jp) - [webserver] - Algorithm based for on-Target sgRNA design.
- [CRISPR-DO](http://cistrome.org/crispr/) - [webserver] - A web application for the design and optimization of guide sequences that target both coding and non-coding regions in spCas9 CRISPR system across human, mouse, zebrafish, fly and worm genomes.
- [CRISPR-FOCUS](http://cistrome.org/crispr-focus/) - [webserver] -  A web-based platform to search and prioritize sgRNAs for CRISPR screening experiments. 
- [CRISPR-P](http://crispr.hzau.edu.cn/CRISPR2/) - [webserver] - One of the most popular tools for sgRNA design in plants with minimal off-target effects.
- [CRISPRscan](http://www.crisprscan.org) - [webserver] - A novel algorithm to predict gRNA efficiency.
- [CRISPRTarget](http://bioanalysis.otago.ac.nz/CRISPRTarget/crispr_analysis.html) - [webserver] - A tool to explore the targets of CRISPR RNAs.
- [CROP-IT](http://cheetah.bioch.virginia.edu.proxygw.wrlc.org/AdliLab/CROP-IT/homepage.html) - [webserver] - A web tool developed to assist biologists in designing CRISPR/CAS9 sgRNAs, by predicting the off-target sites and ranking them according to the chances of occurrence.
- [E-CRISP:](http://www.e-crisp.org/E-CRISP/) - [webserver] - A algorithm is available for twelve organisms and can be easily extended to design both sgRNA and pgRNA.
- [flyCRISPR](http://targetfinder.flycrispr.neuro.brown.edu) - [webserver] - Specificity for Drosophila to find CRISPR target sites and evaluate each identified CRISPR target.
- [GPP Web Portal](https://portals.broadinstitute.org/gpp/public/) - [webserver] -  A web-based platform to generate matching sgRNA knockout (CRISPRko) designs for a mouse or human gene, transcript or target sequence.
- [pgRNADesign](https://bitbucket.org/liulab/pgrnadesign.git) - [python] -  A algorithm to design paired gRNAs for knocking out long non-coding RNAs (lncRNAs).
- [Protospacer Workbench](http://www.protospacer.com) - [Software] - Protospacer Workbench offers an interface for finding, evaluating, and sharing Cas9 guide-RNA (gRNA) designs. 
- [SSC](http://cistrome.org/SSC/) - [webserver] - A sequence model for predicting sgRNA efficiency in CRISPR/Cas9 knockout experiments.
- [WGE](https://www.sanger.ac.uk/htgt/wge/) - [webserver] - A algorithm shows CRISPR sites (paired or single) in and around genes and scores the pairs for potential off-target sites, and browse individual and paired CRISPR sites across human, mouse.
- [WU-CRISPR](http://crispr.wustl.edu) - [webserver] - A web tool to design gRNA for CRISPR/Cas9 Knockout system.


## Genome editing outcomes and predictions

- [CRISPResso](https://github.com/lucapinello/CRISPResso) - [python, webserver] - A software pipeline for the analysis of targeted CRISPR-Cas9 sequencing data. This algorithm allows for the quantification of both non-homologous end joining (NHEJ) and homologous directed repair (HDR) occurrences.
- [CRISPR-GA](http://crispr-ga.net) - [webserver] -  A platform to assess the quality of gene editing using next gen sequencing data to quantify and characterize insertions, deletions and homologous recombination.
- [Microhomology-Predictor](http://www.rgenome.net/mich-calculator/) - [webserver] - A web tool can simply predict the mutation patterns caused by microhomology-mediated end joining (MMEJ) and estimate how frequently unwanted in-frame deletions would be happened.

## Screening analysis

- [BAGEL](https://sourceforge.net/projects/bagel-for-knockout-screens/) - [python] - A algorithm is designed to identify negatively selected genes, by calculating a Bayes factor for each gene representing a confidence measure that the gene knockout results in a fitness defect. Bayesian analysis of gene knockout screens using pooled library CRISPR or RNAi.
- [casTLE](https://bitbucket.org/dmorgens/castle) - [python] - Based on empirical Bayesian framework to account for multiple sources of variability, including reagent efficacy and off-target effects.
- [CaRpools](https://github.com/boutroslab/caRpools) - [R] - A pipeline for end-to-end analysis of pooled CRISPR/Cas9 screening data. Including in-depth analysis of screening quality and sgRNA phenotypes.
- [CERES](https://depmap.org/ceres/) - [R] -  A algorithm to estimate gene-dependency levels from CRISPR-Cas9 essentiality screens while accounting for this effect.
- [CrispRVariants](https://github.com/markrobinsonuzh/CrispRVariants) - [R] - An R-based toolkit for counting, localising and plotting targeted insertion and deletion variant combinations from CRISPR-Cas9 mutagenesis experiments.
- [JACKS](https://github.com/felicityallen/JACKS) [python] - A Bayesian method that jointly analyses screens performed with the same guide RNA library.
- [MAGeCK](https://bitbucket.org/liulab/mageck) - [python] - Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout (MAGeCK) for prioritizing single-guide RNAs, genes and pathways in genome-scale CRISPR/Cas9 knockout screens. 
- [MAGeCK-VISPR](https://bitbucket.org/liulab/mageck-vispr) - [python] - A comprehensive quality control, analysis and visualization workflow for CRISPR/Cas9 screens.
- [MAGeCKFlute](https://bitbucket.org/liulab/mageckflute/) - [R] - A pipeline for perform computational analysis of CRISPR screens. MAGeCKFlute combines the MAGeCK and MAGeCK-VISPR algorithms and incorporates additional downstream analysis functionalities.
- [STARS](https://portals.broadinstitute.org/gpp/public/software/stars) - [python] - A gene-ranking algorithm for genetic perturbation screens, computing a score for genes using the probability mass function of a binomial distribution. To analyze either shRNA or sgRNA based screening data.

