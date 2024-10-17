# Sweden surpasses the UNAIDS 95-95-95 Target: Estimating HIV-1 incidence between 2003-2022

The following codes have been used to generate the results found in the paper "Sweden surpasses the UNAIDS 95-95-95 Target: Estimating HIV-1 incidence between 2003-2022" published in Eurosurveillance on the 10-17-2024.

https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2024.29.42.2400058 (https://doi.org/10.2807/1560-7917.ES.2024.29.42.2400058)

The codes decribed below and available in this repositry require the package **biophybreak** available throught the following link https://github.com/MolEvolEpid/biophybreak, and can be installed into your R libraries repositry.

**prep.incidence.data.new.R** : Takes the required files containing demographic data and biomarker data to construct the dataframe of the required strcutre in order to create the file required to be passed in the incidence estimation.

**incidence.mc.functions.R** : Contains all the required functions to be using in **incidence.estimation.mc.new.R** in order to estimate incidence.

**incidence.estimation.mc.new.R** : Run the incidence model to produce a file that contains the incidence estimations for 10,000 Monte Carolo samples.

**incidence.parse.output.new.gamma.R** : Parses the outputs from **incidence.estimation.mc.new.R** in order to produce plots and mean estimates and confidence intervals.

**incidence.stat.plots.R** : Produces all figures relating to the incidence estimation.

Please note these files will not be maintained by the original author, and may see no maintance at all in the future.
