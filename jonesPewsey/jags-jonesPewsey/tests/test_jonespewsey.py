#!/usr/bin/env python3
"""Tests for the JAGS Jones-Pewsey module.

Checks random sampling, density normalisation, the von Mises limit,
and that psi can be recovered as a continuous parameter.
"""

from __future__ import annotations

import math
import os
import statistics
import subprocess
import tempfile
from pathlib import Path

MODULE_DIR = Path(__file__).resolve().parent.parent
JAGS_LIBS = str(MODULE_DIR)
TWO_PI = 2.0 * math.pi


def run_jags(
    model: str,
    data_lines: list[str],
    *,
    n_iter: int = 5000,
    burnin: int = 500,
    inits_lines: list[str] | None = None,
    monitor: str = "mu",
    stem: str = "out",
) -> list[float]:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        model_file = tmp_path / "model.bug"
        data_file = tmp_path / "data.txt"
        script_file = tmp_path / "script.jag"
        model_file.write_text(model)
        data_file.write_text("\n".join(data_lines) + "\n")

        init_block = ""
        if inits_lines:
            inits_file = tmp_path / "inits.txt"
            inits_file.write_text("\n".join(inits_lines) + "\n")
            init_block = f'parameters in "{inits_file}"\n'

        script = f"""load jonespewsey
model in "{model_file}"
data in "{data_file}"
compile
{init_block}initialize
update {burnin}
monitor {monitor}
update {n_iter}
coda {monitor}, stem("{tmp_path / stem}")
"""
        script_file.write_text(script)
        env = os.environ.copy()
        env["JAGS_LIBS"] = JAGS_LIBS
        result = subprocess.run(
            ["jags", str(script_file)],
            capture_output=True,
            text=True,
            env=env,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                "JAGS failed\n"
                f"stdout:\n{result.stdout}\n"
                f"stderr:\n{result.stderr}"
            )

        chain_file = tmp_path / f"{stem}chain1.txt"
        samples = []
        for line in chain_file.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                try:
                    samples.append(float(parts[1]))
                except ValueError:
                    pass
        return samples


def circular_mean(x: list[float]) -> float:
    s = sum(math.sin(v) for v in x) / len(x)
    c = sum(math.cos(v) for v in x) / len(x)
    return math.atan2(s, c)


def quantile(x: list[float], q: float) -> float:
    ys = sorted(x)
    if not ys:
        raise ValueError("empty sample")
    pos = q * (len(ys) - 1)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return ys[lo]
    w = pos - lo
    return ys[lo] * (1.0 - w) + ys[hi] * w


def linspace(start: float, stop: float, n: int) -> list[float]:
    if n < 2:
        return [start]
    step = (stop - start) / n
    return [start + i * step for i in range(n)]


def run_jags_density(x: list[float], mu: float, kappa: float, psi: float) -> list[float]:
    data_lines = [
        f'"M" <- {len(x)}',
        '"x" <- c(' + ", ".join(f"{v:.16g}" for v in x) + ")",
        f'"mu" <- {mu:.16g}',
        f'"kappa" <- {kappa:.16g}',
        f'"psi" <- {psi:.16g}',
    ]
    return run_jags(
        """model {
  dummy ~ dnorm(0, 1)
  for (i in 1:M) {
    dens[i] <- dJonesPewsey(x[i], mu, kappa, psi)
  }
}""",
        data_lines,
        n_iter=1,
        burnin=0,
        monitor="dens",
        stem="dens",
    )


def main() -> None:
    so = MODULE_DIR / "jonespewsey.so"
    if not so.exists():
        raise SystemExit(f"Build the module first: make -C {MODULE_DIR}")

    print("=== 1. randomSample support and von Mises limit (psi = 0) ===")
    samples = run_jags(
        """model {
  for (i in 1:N) {
    x[i] ~ dJonesPewsey(0, 4, 0)
  }
}""",
        ['"N" <- 4000'],
        n_iter=1,
        burnin=0,
        monitor="x",
        stem="prior0",
    )
    invalid = sum(1 for v in samples if v < 0.0 or v >= TWO_PI)
    print(f"  invalid samples: {invalid} / {len(samples)}")
    assert invalid == 0
    cm = circular_mean(samples)
    print(f"  circular mean: {cm:.3f} (expect ~0)")
    assert abs(cm) < 0.15
    print("  OK")

    print("\n=== 2. density integrates to 1 for continuous psi values ===")
    ngrid = 721
    grid = linspace(0.0, TWO_PI, ngrid)
    dtheta = TWO_PI / ngrid
    for psi in (-1.0, -0.4, 0.0, 0.3, 1.0):
        dens = run_jags_density(grid, mu=0.7, kappa=2.5, psi=psi)
        mass = sum(dens) * dtheta
        print(f"  psi={psi:5.1f}: integral={mass:.5f}")
        assert abs(mass - 1.0) < 0.02, f"density not normalized at psi={psi}"
    print("  OK")

    print("\n=== 3. recover continuous psi from data ===")
    true_psi = -0.55
    x = run_jags(
        """model {
  for (i in 1:N) {
    x[i] ~ dJonesPewsey(mu, kappa, psi)
  }
  mu <- 0.4
  kappa <- 3.0
  psi <- -0.55
}""",
        ['"N" <- 250'],
        n_iter=1,
        burnin=0,
        monitor="x",
        stem="sim",
    )
    x = [math.fmod(v, TWO_PI) for v in x[:250]]
    x = [v + TWO_PI if v < 0 else v for v in x]
    data_lines = [
        f'"N" <- {len(x)}',
        '"x" <- c(' + ", ".join(f"{v:.16g}" for v in x) + ")",
    ]
    psi_post = run_jags(
        """model {
  for (i in 1:N) {
    x[i] ~ dJonesPewsey(mu, kappa, psi)
  }
  mu ~ dnorm(0, 0.001)
  kappa ~ dgamma(0.1, 0.1)
  psi ~ dunif(-2, 2)
}""",
        data_lines,
        inits_lines=[
            f'"mu" <- {circular_mean(x):.16g}',
            '"kappa" <- 2.0',
            '"psi" <- 0.0',
        ],
        n_iter=4000,
        burnin=1000,
        monitor="psi",
        stem="psi",
    )
    psi_mean = statistics.fmean(psi_post)
    psi_lo = quantile(psi_post, 0.05)
    psi_hi = quantile(psi_post, 0.95)
    print(
        f"  true psi={true_psi:.2f}; posterior mean={psi_mean:.2f} "
        f"90% CI=[{psi_lo:.2f}, {psi_hi:.2f}]"
    )
    assert psi_lo < true_psi < psi_hi, "continuous psi not recovered"
    n_unique = len({round(v, 6) for v in psi_post})
    print(f"  unique psi samples (6 d.p.): {n_unique}")
    assert n_unique > 100, "psi looks discrete / not continuously sampled"
    print("  OK")

    print("\n=== 4. lowercase alias djonespewsey ===")
    alias = run_jags(
        """model {
  for (i in 1:N) {
    x[i] ~ djonespewsey(0, 1, 0.5)
  }
}""",
        ['"N" <- 50'],
        n_iter=1,
        burnin=0,
        monitor="x",
        stem="alias",
    )
    assert alias and all(0.0 <= v < TWO_PI for v in alias)
    print("  OK")

    print("\n=== 5. two (kappa, psi) pairs (shared-mu / two set sizes) ===")
    ngrid = 361
    grid = linspace(0.0, TWO_PI, ngrid)
    dtheta = TWO_PI / ngrid
    data_lines = [
        f'"M" <- {ngrid}',
        '"x" <- c(' + ", ".join(f"{v:.16g}" for v in grid) + ")",
    ]
    dens = run_jags(
        """model {
  dummy ~ dnorm(0, 1)
  kappa3 <- 4
  psi3 <- -0.4
  kappa6 <- 1.5
  psi6 <- 0.3
  for (i in 1:M) {
    d3[i] <- dJonesPewsey(x[i], 0, kappa3, psi3)
    d6[i] <- dJonesPewsey(x[i], 1, kappa6, psi6)
  }
}""",
        data_lines,
        n_iter=1,
        burnin=0,
        monitor="d3",
        stem="two3",
    )
    dens6 = run_jags(
        """model {
  dummy ~ dnorm(0, 1)
  kappa3 <- 4
  psi3 <- -0.4
  kappa6 <- 1.5
  psi6 <- 0.3
  for (i in 1:M) {
    d3[i] <- dJonesPewsey(x[i], 0, kappa3, psi3)
    d6[i] <- dJonesPewsey(x[i], 1, kappa6, psi6)
  }
}""",
        data_lines,
        n_iter=1,
        burnin=0,
        monitor="d6",
        stem="two6",
    )
    mass3 = sum(dens) * dtheta
    mass6 = sum(dens6) * dtheta
    print(f"  set-size-3-like integral={mass3:.5f}")
    print(f"  set-size-6-like integral={mass6:.5f}")
    assert abs(mass3 - 1.0) < 0.02
    assert abs(mass6 - 1.0) < 0.02
    print("  OK")

    print("\n=== 6. moderate-large kappa still normalizes ===")
    ngrid_k = 2001
    grid_k = linspace(0.0, TWO_PI, ngrid_k)
    dtheta_k = TWO_PI / ngrid_k
    dens_k = run_jags_density(grid_k, mu=0.0, kappa=8.0, psi=-0.3)
    mass_k = sum(dens_k) * dtheta_k
    print(f"  kappa=8, psi=-0.3: integral={mass_k:.5f}")
    assert abs(mass_k - 1.0) < 0.03
    dens_vm = run_jags_density(grid_k, mu=0.0, kappa=20.0, psi=0.0)
    mass_vm = sum(dens_vm) * dtheta_k
    print(f"  kappa=20, psi=0 (von Mises limit): integral={mass_vm:.5f}")
    assert abs(mass_vm - 1.0) < 0.03
    print("  OK")

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
