### Cluster representation model

Uses the perceptual reproduction representation as a reference and infers the differences to the memory reproduction and similarity comparison representation. The difference weights \gammaPM and \gammaPS measure the overall magnitude of these differences, and can be compared to their prior \gammaPrior to provide Bayes factors using the Savage-Dickey method.

`clusterRepresentation_jags.txt` implements the model as a graphical model in JAGS.

`clusterRepresentation.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is the one used in the paper.