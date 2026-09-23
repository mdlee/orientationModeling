# JAGS Jones-Pewsey module

JAGS extension providing the Jones–Pewsey family of symmetric circular
distributions (Jones & Pewsey, 2005). The shape parameter `psi` is a
**continuous** real parameter, not a discrete grid.

## Distribution

```bugs
x ~ dJonesPewsey(mu, kappa, psi)
```

The lowercase alias `djonespewsey` is also registered.

| Parameter | Description | Support |
| --------- | ----------- | ------- |
| `mu` | Mean direction (radians); any real value, wrapped internally | `(-Inf, Inf)` |
| `kappa` | Concentration | `[0, Inf)` |
| `psi` | Shape (continuous) | `(-Inf, Inf)` |

Special cases of `psi`:

- `psi -> 0`: von Mises
- `psi = -1`: wrapped Cauchy
- `psi = 1`: cardioid
- `kappa = 0` (any `psi`): uniform on `[0, 2*pi)`

Density:

```
f(x | mu, kappa, psi) = (cosh(kappa*psi) + sinh(kappa*psi)*cos(x-mu))^(1/psi)
                        / (2*pi * P_{1/psi}(cosh(kappa*psi)))
```

where `P_nu` is the Legendre function of the first kind. The `psi -> 0`
limit is `exp(kappa * cos(x-mu)) / (2*pi * I_0(kappa))`.

The normalising constant is evaluated numerically (periodic trapezoid
rule) and is a smooth function of continuous `psi`, so `psi` can have
an ordinary JAGS prior, e.g. `psi ~ dunif(-2, 2)` or `psi ~ dnorm(0, 1)`.

## Density function

The same name is also available as a BUGS function (R-style):

```bugs
dens[i] <- dJonesPewsey(x[i], mu, kappa, psi)
```

JAGS also provides `logdensity.JonesPewsey(x, mu, kappa, psi)`.

## Build and install

Requires JAGS development headers and libraries.

```bash
cd jonesPewsey/jags-jonesPewsey
make
sudo make install
```

Run checks with `make test` (uses the local `jonespewsey.so` via `JAGS_LIBS`).

The module is installed to `$(libdir)/JAGS/modules-4/jonespewsey.so`.

The log-normalizer is cached over recent `(kappa, psi)` pairs, so a model with two set-size-specific shapes (shared `mu`) does not recompute the constant on every trial.

## Usage

```r
library(rjags)

model <- "
model {
  for (i in 1:N) {
    theta[i] ~ dJonesPewsey(mu, kappa, psi)
  }
  mu ~ dnorm(0, 0.001)
  kappa ~ dgamma(0.1, 0.1)
  psi ~ dunif(-2, 2)
}
"

jags.model(model, data = list(N = 10), modules = "jonespewsey")
```

From the JAGS command line:

```
load jonespewsey
```

To load from a build directory without installing:

```bash
JAGS_LIBS=/path/to/jags-jonesPewsey jags myscript.jag
```

## Example model

See `examples/jonespewsey.bug`.

## Circular mean estimation

If you put a uniform prior on `mu` over `[0, 2*pi)`, JAGS scalar samplers
treat that interval as a line with hard endpoints. Prefer a wide unbounded
prior such as `mu ~ dnorm(0, 0.001)`, initialize `mu` at the sample circular
mean, and summarize with the circular mean of posterior samples.

## References

- Jones, M. C., & Pewsey, A. (2005). A family of symmetric distributions
  on the circle. *Journal of the American Statistical Association*,
  100(472), 1422–1428.
- Wabersich, D., & Vandekerckhove, J. (2014). Extending JAGS: A tutorial
  on adding custom distributions. *Behavior Research Methods*, 46, 15–28.
