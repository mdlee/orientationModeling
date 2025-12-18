### Memory reproduction model

On each trial 3 or 6 oriented line stimuli are presented then removed. One is then probed as the target, and the participant uses an interface to reproduce it.

`memoryReproduction_jags.txt` and `memoryReproduction_stan.txt` implement the full memory reproduction comparison model as a graphical model in JAGS and Stan, respectively.

`memoryReproductionNoSwap_jags.txt` and `memoryReproductionNoSwap_stan.txt` implement  a reduced version of the model that removes the swap error process.

`memoryReproduction.m` and `memoryReproductionNoSwap.m` are MATLAB scripts that apply the models to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.

`memoryReproduction.Rmd` is a RMarkdown script that applies the memory reproduction with swap errors model to data using R2jags. Instructions on how to run the script are provided within the file. This file can be easily changed to run the alternative memory reproduction model (no swap errors).

`memoryReproduction.html` is the HTML output of the RMarkdown script. This can be opened within any browser to review the results of the models.