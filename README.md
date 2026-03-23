# clonal_structure

Python code, experimental and derived datasets for "Universal bacterial clade dynamics dominate under predation despite altered phenotypes and mutation targets" study (now published in Evolution; DOI: https://doi.org/10.1093/evolut/qpag043).

The description below is adapted from the Dryad archive for the corresponding study (Dataset DOI: 10.5061/dryad.ffbg79d8s).

Description of the data and file structure
Experimental and derived datasets for "Universal bacterial clade dynamics dominate under predation despite altered phenotypes and mutation targets" study. Experimental data provided includes LogPhase measurements of each of the five evolved-alone/coevolved-with-predator single species experiments, growing on the three different growth media (files in "phenotyping" directory; see full article for details). Genomic variants are provided in .vcf files ("vcfs" directory). Raw genomic sequences are deposited in ENA (accession: PRJEB85532). The remaining files are program code and various intermediate data to facilitate result replicability.

Files and variables
File: chosen_sample_df.csv
Description: a table containing a list of curated genomic sample filenames for each experimental population.

File: clonal_order_df.csv
Description: a table containing a list of tuples for each experimental population. Each tuple contains the number of emerged clone in the population and its emergence time.

File: cluster_clone_match.tsv
Description: a table containing clone number, clone name (represented in both a single character or a string of characters, following clone hierarchy), and the time of its emergence.

File: rep_time.tsv
Description: a table containing a list for each experimental population, representing the curated list of sampling time points in days during the long-term experimental evolution experiment.

File: Sampling_experiment_day_matches.txt
Description: the main experimental metadata file.

File: sequencing_sample_information.csv
Description: metadata for sequencing samples, includes a column denoting if a sample was selected for clonal structure inference.

Variables
id: unique experimental population identifier;
bacterial_strain: species name;
organism_type: bacteria or predator;
experiment: evolutionary history identifier, denoting evolving-alone (NP) and coevolving with predator (PS) populations;
replicate: replicate experimental population (three for each species and evolutionary history combination);
sampling_date: the calendar date the experimental population was sampled;
day_in_experiment: the day in experiment that corresponds to the sampling date;
unreliable_pop_size: a boolean variable denoting whether the experimental population size (prey optical density) was trustworthy (e.g., the optical density close to the detection threshold was considered unreliable);
notes: additional notes related to the experimental population;
OD600: optical density measurement of the prey population in the experimental population;
pred_cells_ml: ciliate cell counts in the experimental population;
pred_corrected: corrected cilliate cell counts in the experimental population.
File: clonal_structure-github.zip
Description: a cloned GitHub repository (contains the same code and data as presented here).

File: pyclone_inputs.tar.gz
Description: an archive with input files to PyClone software for each experimental population, representing curated genomic variant trajectories.

File: pyclone_output.tar.gz
Description: an archive with PyClone software output files for each experimental population.

File: phenotyping.tar.gz
Description: an archive with LogPhase measurement device output files, measuring optical density of the prey populations, sampled at different time points across the long-term experimental evolution experiment, growing on three different media (see main article for details).

File: vcfs.tar.gz
Description: an archive containing the detected genomic variants (vcf files) for each experimental population.

File: ancestral_evolved_aucs.csv
Description: a table with computed ancestral to (co-)evolved prey growth area under the curve (AUC) ratios, used to compute linear mixed models in R.

File: no_interaction_R2_decomposition.csv
Description: linear mixed model R squared decomposition results from R, to be used with Python code.

Code/software
All code written in Python (version 3.10.14) and R (version 4.5.2). This repository includes program code files listed below.

Jupyter iPython notebooks:

clonal_structure_inference.ipynb - this notebook shows how the clonal structure is computed for a given experiment.

mutation_statistics.ipynb - here we compute mutation count and recurrence statistics and make figures reported in the main article.

phenotype_analysis.ipynb - here we analyse our phenotypic data and make figures reported in the main article.

tree_analysis.ipynb - here we analyse clonal structure trees, compute statistics and make figures reported in the main article.

Python script files:

clonal_structure.py - main code with clonal structure inference algorithm.

utils.py - some helper functions.

R markdown files:

linear_mixed_models.Rmd - here we utilise linear mixed models to analyse growth phenotype (ancestral to (co-)evolved prey growth AUC ratios) measurements of our experimental species growing on three different media.

Access information
Other publicly accessible locations of the data:

Dryad: 10.5061/dryad.ffbg79d8s
Data was derived from the following sources:

ENA (accession: PRJEB85532)
