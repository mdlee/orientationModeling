### Similarity comparison model

On each trial, two pairs of oriented line stimuli are presented. The response is a binary decision about which of the two pairs has more similar stimuli.

`similarityComparison_jags.txt` implements the similarity comparison model as a graphical model in JAGS. A direct Stan implementation is difficult, if not impossible, because of the use of censoring.

`similarityComparison.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.

`similarityComparison.Rmd` is a RMarkdown script that applies the model to data using R2jags. Instructions on how to run the script are provided within the file.

`similarityComparison.html` is the HTML output of the RMarkdown script. This can be opened within any browser to review the results of the models.