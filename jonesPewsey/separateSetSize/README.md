# Jones–Pewsey memory fits, one analysis per set size

These scripts fit the memory-reproduction Jones–Pewsey models **separately** for set size 3 and set size 6. Each run has its own psychological map `mu`, concentration `sigma` / `kappa`, shape `psi`, and (for the swap model) swap weights `omega`.

That is a different model from the paper-style specification in the parent folder, which uses **one shared `mu`** for both set sizes.

## When to use this folder

Use these fits when you want set-size-specific representations, or as a simpler sampler (one `(kappa, psi)` pair per JAGS run, like `perceptualReproductionJP`).

Use the parent-folder JP drivers when you want the original joint model:

| Shared-`mu` model (parent folder) | This folder |
| -------------------------------- | ----------- |
| `../memoryReproductionNoSwapJP.m` | `noSwap.m` |
| `../memoryReproductionJP.m` | `withSwap.m` |

## What is independent, and what is not

In the joint model, **noise and swap parameters are already set-size-specific** (`sigma`, `psi`, `omega`). The only node that couples the two set sizes is `mu`. Splitting the data therefore:

- leaves `sigma` / `psi` / `omega` with the same meaning
- **stops pooling** the stimulus map across set sizes (each set size gets its own `mu`)

## Scripts

| File | Model | JAGS |
| ---- | ----- | ---- |
| `noSwap.m` | Target always recalled | `noSwap_jags.txt` |
| `withSwap.m` | Mixture over similarity-ranked items | `withSwap_jags.txt` |

Each driver loops over set sizes 3 then 6.

## Outputs

Written under this folder:

- `storage/noSwap_setSize3_tomicBays_jags.mat`
- `storage/noSwap_setSize6_tomicBays_jags.mat`
- `storage/withSwap_setSize3_tomicBays_jags.mat`
- `storage/withSwap_setSize6_tomicBays_jags.mat`
- `figures/tomicBays_noSwap_setSize3.png` (and `.eps`, plus set size 6 / withSwap analogues)

## Requirements

Same as the other JP drivers: build `../jags-jonesPewsey` (`make`), MATLAB with trinity / `orientationMcmcProtocol`. Run the `.m` files from this directory (they `cd` to themselves).
