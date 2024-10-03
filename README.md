# Multicountry genomic analysis underscores regional cholera spread in Africa
**Authors**: Gerald Mboowa, Nathaniel Lucero Matteson, Collins Kipngetich Tanui, Mpanga Kasonde, Guyguy Kusanzangana Kamwiziku, Olusola Anuoluwapo Akanbi, Jucunú Johane Elias Chitio, Mathews Kagoli, Rene Ghislain Essomba, Alisen Ayitewala, Isaac Ssewanyana, Blaise Mboringong Akenji, Adrienne Aziza Amuri, Andrew S Azman, Olajumoke Atinuke Babatunde, Valentina Josiane Ngo Bitoungui, Espoir Malembaka Bwenge, Francis Ongole, Chimaobi Emmanuel Chukwu, Nália Ismael, Otridah Kapona, Osvaldo Laurindo, Placide Kingebeni Mbala, Georges Alain Etoundi Mballa, Imelda Carlos Zulfa Miambo, Alex Ansaye Mwanyongo, Grace Najjuka, Joseph Mutale, Kunda Musonda, Allan Muruta Niyonzima, Mirriam Ethel Nyenje, Michael Popoola, Doreen Shempela, Christiane Medi Sike, Sofião Manjor Sitoe, Dorcas Waruguru Wanjohi, Placide Okitayemba Welo, Mtisunge Yelewa, Sebastian Yennan, Lucius Ziba, CholGen Consortium, Joseph Ephram Bitilinyu-Bangoh, Roma Chilengi, Hamsatou Hadja, Jide Idris, José Paulo Maurício Langa, Daniel Bamuleka Mukadi, Susan Nabadda, Amanda K Debes, David A Sack, Jean Kaseya, Yenew Kebede Tebeje, Shirlee Wohl, Sofonias Kifle Tessema

DOI: TBD
## Abstract
Cholera remains a significant public health burden in many countries in sub-Saharan Africa, though the exact mechanisms of bacterial emergence and spread remain undefined. 
We used genomic data from 768 Vibrio cholerae O1 isolates predominantly collected between 2019-2024 to generate the largest dataset of V. cholerae genomes sequenced locally in Africa, which we used to interrogate recent patterns of spread, including the rapid circulation of the AFR15 lineage associated with unusually large outbreaks in Southern Africa. 
We provide evidence for the movement of this lineage into new African Member States and confirm differences in transmission dynamics previously observed in West versus East Africa, though cross-border transmission is prevalent in both regions of the continent. 
Despite observed differences, evolutionary processes are similar across lineages and we find no evidence for significant changes in antimicrobial resistance genotypes. 
Taken together, our findings point to the importance of regionally coordinated cross-border surveillance and implementation of interventions across country borders, and highlight the value of locally generated genomic data for answering crucial questions about the spread of cholera in Africa.

## Organization of this repository
This repository contains data and code necessary to reproduce the results in our paper. 
- The `data/` directory contains the sequences and metadata generated by this study as well as the background dataset used for the phylogenetic analysis. Shape files for visualization are also present here.
- The `analysis/` directory contains the code for generating the figures of our paper. Most figures/panels have their own notebook, though occassionally they are merged into one. Child directorys contain the code necessary for transforming the data (for example, extracting Markov jumps). Individual figure panels are saved in the `plots/` subdirectory.
- The `beast-analyses/` directory contains the XMLs, log files, and posterior trees from our BEAST analyses.
- The `figures/` directory contains the final full resolution figures for our paper.
- The `nextstrain/` directory contains the nextstrain build of the third wave 7PET phylogeny we analyzed as part of this manuscript. You can access the build [on Nextstrain](https://nextstrain.org/community/CholGen/Regional-Analysis-2024_wave3) or with a local auspice server. 

## Reproducing results
Most the results of our paper can be reproduced using the included Snakemake pipeline. 
The pipeline will re-execute all the notebooks in the `analysis/` directory. 
You can install the prerequisites for our analyses using the provided [conda environment](environment.yaml) (though we recommend using the [mamba wrapper](https://github.com/mamba-org/mamba)).
Installing the prerequisities should take less than 10 minutes.
Once the environment is installed and activate, the full pipeline can be run with 
```bash
snakemake --force-all --rerun-incomplete --cores <number-of-cores>
```
The snakemake pipeline should take less than 30 minutes to regenerate all figures.

The only results that will not be reproduced by this script are the BEAST runs.
However, these can be rerun with the provided XMLs, but note that these runs take a considerable amount of time to complete (>300 hours).

