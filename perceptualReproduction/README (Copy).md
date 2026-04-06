### Perceptual reproduction model

On each trial an oriented line stimulus is presented and the participant uses an interface to reproduce it.

`perceptualReproduction_jags.txt` implement the perceptual reproduction comparison model as a graphical model in JAGS.

<<<<<<< HEAD
`perceptualReproduction.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented in the paper.

The JAGS script implements the circular normal distribution using a latent variable that allows any stimulus location to be shifted by +pi or -pi. We attempted to use the von Mises module [here](https://github.com/yeagle/jags-vonmises) for JAGS, but believe it does not work well. It did not produce sensible results even for very simple inference problems. We would be pleased to hear from anybody who is able to provide a working von Mises implementation of the model in JAGS.

Given difficulties in producing a Stan implementation of the memory reproduction and similarity comparisons (as documented in their READMEs), we did not pursue a Stan implementation of the perceptual reproduction model, although we think it would be straightforward using the inbuilt von Mises distribution.
=======
`perceptualReproduction.m` is a MATLAB script that applies the model to data. The `tomicBays` empirical data set from [here](https://psycnet.apa.org/record/2023-21056-001) is implemented, but others could be added.

`perceptualReproduction.Rmd` is a RMarkdown script that applies the model to data using R2jags. Instructions on how to run the script are provided within the file.

`perceptualReproduction.html` is the HTML output of the RMarkdown script. This can be opened within any browser to review the results of the models.
>>>>>>> 77e93d658406209d1d011dd70d30b472b697abc0
