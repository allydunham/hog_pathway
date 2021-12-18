# *S. cerevisiae* phenotype prediction analysis.

Scripts analysing the growth of 1,011 *S. cerevisiae* strains and attempting to predict growth phenotypes from genotype.
Growth phenotypes were sourced from [Peter et al. (2018)](https://www.nature.com/articles/s41586-018-0030-5) and measured by the Beltrao lab, where I performed this anlysis as part of my PhD.
The repo also contains follow up analysis on the growth phenotypes of 4 different gene deletion libraries, which is included in [Galardini et al. (2019)](https://onlinelibrary.wiley.com/doi/abs/10.15252/msb.20198831).
I have also performed another genotype to phenotype analysis on the same strains, incorporating transcriptomic and proteomic data ([repo](https://github.com/allydunham/yeast_strains)).
All of these analyses make use of the P(aff) score to model genotype, which is detailed in [Jelier et al. (2011)](https://www.nature.com/articles/ng.1007).


There are five main analyses:

* Link between genetic difference and phenotype difference
* Associations between genes and significant growth phenotype
* Modelling phenotype based on genotype of sets of associated genes
* Probability model of the HOG pathway
* Comparison of different knockout library phenotype, identifying cases where deletion phenotypes differ based on genetic background

Most of the scripts in this repo rely on specific data that is not all easily available, meaning any researcher seeking ot perform a similar analysis will likely find it easier to work on their own data, simply using these analyses as an outline.
The `pathway.py` module is the exception, implementing a general framework for a probability model predicting pathway activity, which can be easily applied to other pathways.

The data processing and analysis scripts, and modules supporting them are all included in the `bin` directory.
The most up-to-date and well formatted analyses are included in the two `thesis_figure` scripts.
`meta` contains information about genes, particularly those making up the HOG pathway.
`mutfunc` includes two scripts used for processing the raw data from the [Mutfunc]() webserver, which are not generally useful because this data is not publically available.
