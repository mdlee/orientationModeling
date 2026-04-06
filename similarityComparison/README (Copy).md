### Similarity comparison model

On each trial, two pairs of oriented line stimuli are presented. The response is a binary decision about which of the two pairs has more similar stimuli.

`similarityComparison_jags.txt` implements the similarity comparison model as a graphical model in JAGS. There is no need for the latent mixture approach to implementing the circular normal, because ...

<<<<<<< HEAD
`similarityComparison.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is the one used in the paper.

We were not able to implement the model in Stan, because of the use of censoring. In an attempt to produce a Stan implementation, we explored an alternative version of the model that did not use censoring, but instead used an error-of-execution (trembling hand) approach to make the deterministic decision rule stochastic. We were again not able to implement this model scalably in Stan, because it insists on tracking the intermediate mental samples.  We would be pleased to hear from anybody who is able to provide a Stan implementation of either version of the model.
=======
`similarityComparison.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.

`similarityComparison.Rmd` is a RMarkdown script that applies the model to data using R2jags. Instructions on how to run the script are provided within the file.

`similarityComparison.html` is the HTML output of the RMarkdown script. This can be opened within any browser to review the results of the models.
>>>>>>> 77e93d658406209d1d011dd70d30b472b697abc0
