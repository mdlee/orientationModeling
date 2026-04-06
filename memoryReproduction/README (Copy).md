### Memory reproduction model(s)

On each trial 3 or 6 oriented line stimuli are presented then removed. One is then probed as the target, and the participant uses an interface to reproduce it.

`memoryReproduction_jags.txt` implements the full memory reproduction comparison model as a graphical model in JAGS.

`memoryReproductionNoSwap_jags.txt` implements a reduced version of the model that removes the swap error process. This is effectively the perceptual reproduction model.

<<<<<<< HEAD
`memoryReproduction.m` and `memoryReproductionNoSwap.m` are MATLAB scripts that apply the models to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is the one used in the paper.

The JAGS script implements the circular normal distribution using a latent variable that allows any stimulus location to be shifted by +pi or -pi. We attempted to use the von Mises module [here](https://github.com/yeagle/jags-vonmises) for JAGS, but believe it does not work well. It did not produce sensible results even for very simple inference problems. We would be pleased to hear from anybody who is able to provide a working von Mises implementation of the model in JAGS.

We also attempted a Stan implementation, use the inbuilt von Mises distribution. This does appear to work, but we could not produce an effective implementation of the memory reproduction model, because of its reliance on the discrete latent mixture variables \xi. Our implementation mixed over these variables in the standard way, but failed to converge even with very long MCMC runs with extensive thinning. We suspect an effective Stan implementation is possible for the perceptual reproduction model, but did not pursue it. We would be pleased to hear from anybody who is able to provide a Stan implementation of the memory reproduction model.
=======
`memoryReproduction.m` and `memoryReproductionNoSwap.m` are MATLAB scripts that apply the models to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.

`memoryReproduction.Rmd` is a RMarkdown script that applies the memory reproduction with swap errors model to data using R2jags. Instructions on how to run the script are provided within the file. This file can be easily changed to run the alternative memory reproduction model (no swap errors).

`memoryReproduction.html` is the HTML output of the RMarkdown script. This can be opened within any browser to review the results of the models.
>>>>>>> 77e93d658406209d1d011dd70d30b472b697abc0
