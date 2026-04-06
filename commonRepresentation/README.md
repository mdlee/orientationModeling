### Common representation model

The model assumes a single latent representation generates the perceptual reproduction, memory reproduction, and similarity comparison data. It uses the same decision processes as for the individual models.

`commonRepresentation_jags.txt` implements the common representation model as a graphical model in JAGS.

`commonRepresentation.m` is a MATLAB script that apply the models to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is the one used in the paper.