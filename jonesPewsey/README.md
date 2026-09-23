# Jones–Pewsey models

Jones–Pewsey JAGS module and task drivers, at the same level as `perceptualReproduction/`, `memoryReproduction/`, and `similarityComparison/`. Wrap-copy circular-normal and von Mises drivers stay in those original folders.

Mental variables are Jones–Pewsey on the doubled circle (`[0, 2π)`), with `yCirc = mod(2*y, 2π)`, `phiMu = 2*mu`, and `kappa = 1/sigma^2`. The shape `psi` is a continuous parameter (not an EM grid). Special cases: `psi → 0` von Mises, `psi = -1` wrapped Cauchy, `psi = 1` cardioid.

## Layout

```
jonesPewsey/
├── jags-jonesPewsey/              # C++ JAGS module (x ~ dJonesPewsey)
├── perceptualReproductionJP.m
├── memoryReproductionJP.m         # shared μ; set-size-specific σ, ψ, ω
├── memoryReproductionNoSwapJP.m
├── similarityComparisonJP.m
├── separateSetSize/               # each set size has its own μ
├── figures/
└── storage/                        # gitignored MCMC output
```

## Build the module

```bash
cd jonesPewsey/jags-jonesPewsey && make
```

MATLAB drivers set `JAGS_LIBS` to `./jags-jonesPewsey` and load `cfg.modules = {'jonespewsey'}`. Details: `jags-jonesPewsey/README.md`.

## Drivers

Run from this folder (each script `cd`s to itself). `preLoad = true` loads cached `storage/*.mat` when present.

| Driver | Model |
|--------|--------|
| `perceptualReproductionJP.m` | One `σ`, one `ψ` |
| `memoryReproductionJP.m` | Shared `μ`; set-size-specific `σ`, `ψ`, `ω` |
| `memoryReproductionNoSwapJP.m` | Shared `μ`; no swap |
| `similarityComparisonJP.m` | One `σ`, `ψ`; four latent JP samples + `dinterval` |
| `separateSetSize/withSwap.m` | Separate `μ` per set size, with swap |
| `separateSetSize/noSwap.m` | Separate `μ` per set size, no swap |

Data still come from `../data`; shared helpers from `../supportingFiles`.
