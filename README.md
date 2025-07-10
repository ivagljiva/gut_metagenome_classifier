# A logistic regression classifier for human gut metagenomes

This repository provides convenient access to the classifier described in the _eLife_ paper "Microbes with higher metabolic independence are enriched in human gut microbiomes under stress" (https://doi.org/10.7554/eLife.89862), which classifies gut metagenome samples as 'IBD' (Inflammatory Bowel Disease, or more generically, 'stressed') vs 'HEALTHY'. 

The classifier works with per-population copy numbers of specific KEGG modules to predict whether a given metagenome sample comes from a healthy person, or a person with a 'stressed gut environment' (from the perspective of the gut microbes) -- in particular, we trained the classifier on samples from individuals with Inflammatory Bowel Disease (IBD), but our tests showed it also worked to identify more generic stress on the gut microbiome in the form of antibiotic treatment. The specific KEGG modules we worked with are those that we found to be enriched in samples from individuals with IBD. See the [paper](https://doi.org/10.7554/eLife.89862) for more details.

The classifier code and necessary input files for training were originally made available as part of the [reproducible workflow](https://merenlab.org/data/ibd-gut-metabolism/) (and its corresponding datapack) for the paper, but to make it more accessible, I have created a script that allows users to run the model directly on their own input files. For clarity, I have also copied the original classifier Jupyter notebook (`metagenome_classifier.ipynb`) and relevant input files here.

If you want to learn how I trained and tested the model, you can run (or just read) the Jupyter notebook. If you want to run the classifier on your own data, read on.

## Installing Dependencies

I recommend using `conda` to establish an environment in which to run the classifier:

```bash
conda create -n gut_classifier python=3.10
conda activate gut_classifier
pip install -r requirements.txt
```

## Checking options

The classifier script accepts various parameters. You can check them by running the following command:

```bash
python run_classifier.py -h
```

## Running the classifier

You need per-population copy numbers (PPCN values) for each of the IBD-enriched modules in each of your metagenome samples. Alternatively, you can work with raw copy numbers (an output of the [anvi'o](https://anvio.org/) program `anvi-estimate-metabolism`) and provide a table of population sizes to normalize with -- here, population size refers to the number of microbial genomes estimated to be present in a given metagenome assembly. There is a list of the IBD-enriched modules that are used as model features [here in this repository](https://github.com/ivagljiva/gut_metagenome_classifier/blob/main/TRAINING_DATA/IBD_ENRICHED_MODULES.txt). Those investigating a different system might require selecting different features and training a new classifier, and in that case I recommend looking at our feature selection, training, and testing strategies in the Jupyter Notebook.

The [reproducible workflow](https://merenlab.org/data/ibd-gut-metabolism/) for the paper describes how we generate these data files, and I recommend following a similar strategy to obtain your input data.

Below, I show basic examples for two different ways to run the classifier script. These examples use data from the `TESTING_DATA` directory, which originate from the antibiotic time-series dataset described in the study, "Recovery of gut microbiota of healthy adults following antibiotic exposure" by Palleja et al (DOI: [10.1038/s41564-018-0257-9](https://doi.org/10.1038/s41564-018-0257-9)). More details can be found in the paper and the reproducible workflow.

Given a tab-delimited table of PPCN values, this is how you run the classifier:

```bash
python run_classifier.py --ppcn-table TESTING_DATA/PALLEJA_ET_AL_PPCN_matrix.txt
```

Given a tab-delimited table of module copy numbers (unnormalized) and a table of per-sample population sizes, this is how you run the classifier:

```bash
python run_classifier.py --copy-numbers TESTING_DATA/antibiotics-module_stepwise_copy_number-MATRIX.txt --populations TESTING_DATA/02_PALLEJA_SAMPLES_INFO.txt
```

Check the program help output for more options.
