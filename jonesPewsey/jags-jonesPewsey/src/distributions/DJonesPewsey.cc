#include "DJonesPewsey.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <rng/RNG.h>
#include <util/nainf.h>

#include "../JonesPewseyMath.h"

namespace jags {
namespace jonespewsey {

namespace {

double sampleVonMises(double kappa, RNG *rng) {
  if (kappa <= 0.0) {
    return kTwoPi * rng->uniform();
  }

  double theta = 0.0;
  if (kappa <= 0.5) {
    while (true) {
      const double u1 = kTwoPi * rng->uniform();
      const double u2 = rng->uniform();
      if (std::log(u2) <= kappa * (std::cos(u1) - 1.0)) {
        theta = u1;
        break;
      }
    }
  } else {
    const double tau = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    const double rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * kappa);
    const double r = (1.0 + rho * rho) / (2.0 * rho);

    while (true) {
      const double z = std::cos(M_PI * rng->uniform());
      const double f = (1.0 + r * z) / (r + z);
      const double c = kappa * (r - f);
      const double u2 = rng->uniform();

      if ((c * (2.0 - c) - u2) > 0.0 ||
          (std::log(c / u2) + 1.0 - c) > 0.0) {
        const double u3 = rng->uniform();
        theta = (std::floor(u3 + 0.5) * 2.0 - 1.0) * std::acos(f);
        break;
      }
    }
  }
  return wrapAngle(theta);
}

double sampleWrappedCauchy(double rho, RNG *rng) {
  if (rho <= 0.0) {
    return kTwoPi * rng->uniform();
  }
  if (rho >= 1.0) {
    return 0.0;
  }
  const double u = rng->uniform();
  const double scale = (1.0 - rho) / (1.0 + rho);
  return wrapAngle(2.0 * std::atan(scale * std::tan(M_PI * (u - 0.5))));
}

}  // namespace

DJonesPewsey::DJonesPewsey(std::string const &name)
    : ScalarDist(name, 3, DIST_SPECIAL) {}

double DJonesPewsey::logDensity(double x, PDFType type,
                                std::vector<double const *> const &par,
                                double const * /*lower*/,
                                double const * /*upper*/) const {
  const double mu = *par[0];
  const double kappa = *par[1];
  const double psi = *par[2];

  if (kappa < 0.0) {
    return JAGS_NEGINF;
  }

  const double c = std::cos(wrapAngle(x) - wrapAngle(mu));
  const double kernel = logUnnormalized(kappa, psi, c);

  if (type == PDF_PRIOR) {
    return kernel;
  }

  return kernel - cachedLogNormalizer(kappa, psi);
}

double DJonesPewsey::l(std::vector<double const *> const & /*parameters*/) const {
  return 0.0;
}

double DJonesPewsey::u(std::vector<double const *> const & /*parameters*/) const {
  return kTwoPi;
}

double DJonesPewsey::randomSample(std::vector<double const *> const &par,
                                  double const * /*lower*/,
                                  double const * /*upper*/,
                                  RNG *rng) const {
  const double mu = *par[0];
  const double kappa = *par[1];
  const double psi = *par[2];

  if (kappa <= 0.0) {
    return kTwoPi * rng->uniform();
  }

  // psi -> 0: von Mises (Best-Fisher).
  if (std::fabs(psi) < kPsiVonMises) {
    return wrapAngle(sampleVonMises(kappa, rng) + mu);
  }

  // psi = -1: wrapped Cauchy with rho = tanh(kappa / 2).
  if (std::fabs(psi + 1.0) < 1.0e-12) {
    return wrapAngle(sampleWrappedCauchy(std::tanh(0.5 * kappa), rng) + mu);
  }

  // Uniform envelope: unnormalized density is at most exp(kappa).
  // For large kappa use a von Mises proposal with empirical-supremum bound.
  const bool use_vm_proposal = (kappa > 0.5 && std::fabs(psi) < 2.0);

  if (!use_vm_proposal) {
    while (true) {
      const double x = kTwoPi * rng->uniform();
      const double u = rng->uniform();
      const double logh = logUnnormalized(kappa, psi, std::cos(x));
      if (std::log(u) <= logh - kappa) {
        return wrapAngle(x + mu);
      }
    }
  }

  double bound = 1.25;
  while (true) {
    const double x = sampleVonMises(kappa, rng);
    const double c = std::cos(x);
    const double log_ratio =
        logUnnormalized(kappa, psi, c) - kappa * c;
    const double ratio = std::exp(std::min(log_ratio, 700.0));
    if (ratio > bound) {
      bound = ratio * 1.0000001;
      continue;
    }
    if (rng->uniform() * bound <= ratio) {
      return wrapAngle(x + mu);
    }
  }
}

double DJonesPewsey::typicalValue(std::vector<double const *> const &par,
                                  double const * /*lower*/,
                                  double const * /*upper*/) const {
  return wrapAngle(*par[0]);
}

bool DJonesPewsey::checkParameterValue(
    std::vector<double const *> const &par) const {
  const double kappa = *par[1];
  const double psi = *par[2];
  return std::isfinite(*par[0]) && std::isfinite(kappa) &&
         std::isfinite(psi) && kappa >= 0.0;
}

bool DJonesPewsey::isSupportFixed(std::vector<bool> const & /*fixmask*/) const {
  return true;
}

}  // namespace jonespewsey
}  // namespace jags
