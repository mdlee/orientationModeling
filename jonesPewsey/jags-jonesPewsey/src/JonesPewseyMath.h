#ifndef JONESPEWSEY_MATH_H_
#define JONESPEWSEY_MATH_H_

#include <algorithm>
#include <cmath>
#include <limits>

#include <JRmath.h>

namespace jags {
namespace jonespewsey {

const double kTwoPi = 6.283185307179586476925286766559;
const double kLogTwoPi = 1.837877066409345483560659472811;
const double kPsiVonMises = 1.0e-8;

inline double wrapAngle(double x) {
  x = std::fmod(x, kTwoPi);
  if (x < 0.0) {
    x += kTwoPi;
  }
  return x;
}

// log(cosh(a) + sinh(a) * c) with c = cos(phi) in [-1, 1].
inline double logCoshSinhCombo(double a, double c) {
  c = std::max(-1.0, std::min(1.0, c));
  if (a >= 0.0) {
    const double t = std::exp(-2.0 * a);
    const double inner = (1.0 + c) + (1.0 - c) * t;
    if (inner > 0.0) {
      return a - std::log(2.0) + std::log(inner);
    }
    return -a - std::log(2.0) + std::log(std::max(1.0e-300, 1.0 - c));
  }
  const double t = std::exp(2.0 * a);
  const double inner = (1.0 - c) + (1.0 + c) * t;
  if (inner > 0.0) {
    return -a - std::log(2.0) + std::log(inner);
  }
  return a - std::log(2.0) + std::log(std::max(1.0e-300, 1.0 + c));
}

inline double logUnnormalized(double kappa, double psi, double cos_phi) {
  if (kappa <= 0.0) {
    return 0.0;
  }
  cos_phi = std::max(-1.0, std::min(1.0, cos_phi));
  if (std::fabs(psi) < kPsiVonMises) {
    return kappa * cos_phi;
  }
  const double log_base = logCoshSinhCombo(kappa * psi, cos_phi);
  if (!std::isfinite(log_base)) {
    return (psi > 0.0) ? -kappa : kappa;
  }
  const double val = log_base / psi;
  if (!std::isfinite(val)) {
    return (psi > 0.0) ? -kappa : kappa;
  }
  return val;
}

inline double logBesselI0(double kappa) {
  if (kappa <= 0.0) {
    return 0.0;
  }
  const double scaled = bessel_i(kappa, 0.0, 2.0);
  if (scaled > 0.0 && std::isfinite(scaled)) {
    return kappa + std::log(scaled);
  }
  const double unscaled = bessel_i(kappa, 0.0, 1.0);
  if (unscaled > 0.0 && std::isfinite(unscaled)) {
    return std::log(unscaled);
  }
  return kappa - 0.5 * std::log(kTwoPi * kappa);
}

// 2F1(a, b; 1; w) for |w| < 1. Returns NaN if the series fails.
inline double hyp2f1_c1(double a, double b, double w) {
  if (!(std::fabs(w) < 1.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (std::fabs(w) < 1.0e-16) {
    return 1.0;
  }
  double term = 1.0;
  double sum = 1.0;
  for (int n = 0; n < 8000; ++n) {
    const double np1 = static_cast<double>(n + 1);
    term *= (a + n) * (b + n) / (np1 * np1) * w;
    sum += term;
    if (!std::isfinite(sum) || !std::isfinite(term)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::fabs(term) < 1.0e-16 * (1.0 + std::fabs(sum))) {
      return sum;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

inline double logCosh(double a) {
  const double aa = std::fabs(a);
  return aa + std::log((1.0 + std::exp(-2.0 * aa)) / 2.0);
}

// log(2*pi*P_{1/psi}(cosh(kappa*psi))) via
// P_nu(z) = z^nu * 2F1(-nu/2, (1-nu)/2; 1; tanh^2(kappa*psi)), nu = 1/psi.
inline double logNormalizer2F1(double kappa, double psi) {
  if (std::fabs(psi) < 1.0e-4) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double a = kappa * psi;
  const double th = std::tanh(a);
  const double w = th * th;
  if (!(w < 0.99)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double aa = -1.0 / (2.0 * psi);
  const double bb = 0.5 * (1.0 - 1.0 / psi);
  const double f = hyp2f1_c1(aa, bb, w);
  if (!(f > 0.0) || !std::isfinite(f)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return kLogTwoPi + logCosh(a) / psi + std::log(f);
}

// Periodic trapezoid on [-pi, pi] (full circle). Grid size grows with
// sqrt(kappa) so the mode is resolved without a 65k-point sweep.
inline double logNormalizerQuad(double kappa, double psi) {
  int M = 512 + static_cast<int>(80.0 * std::sqrt(1.0 + kappa));
  if (kappa > 4.0) {
    M += static_cast<int>(16.0 * kappa);
  }
  if (M > 8192) {
    M = 8192;
  }
  if (M % 2 != 0) {
    ++M;
  }
  const double dtheta = kTwoPi / static_cast<double>(M);
  double sum = 0.0;
  for (int j = 0; j < M; ++j) {
    const double theta = dtheta * static_cast<double>(j);
    const double g = logUnnormalized(kappa, psi, std::cos(theta)) - kappa;
    if (g > -700.0) {
      sum += std::exp(g);
    }
  }
  if (!(sum > 0.0) || !std::isfinite(sum)) {
    return kLogTwoPi + logBesselI0(kappa);
  }
  return kappa + std::log(sum) + std::log(dtheta);
}

inline double logNormalizer(double kappa, double psi) {
  if (kappa <= 0.0) {
    return kLogTwoPi;
  }
  if (std::fabs(psi) < kPsiVonMises) {
    return kLogTwoPi + logBesselI0(kappa);
  }
  const double z2f1 = logNormalizer2F1(kappa, psi);
  if (std::isfinite(z2f1)) {
    return z2f1;
  }
  return logNormalizerQuad(kappa, psi);
}

// Move-to-front cache so two (or more) (kappa, psi) pairs — e.g. two set
// sizes sharing a likelihood — do not recompute the normalizer every trial.
inline double cachedLogNormalizer(double kappa, double psi) {
  const int kCacheN = 32;
  struct Entry {
    double kappa;
    double psi;
    double logZ;
  };
  thread_local Entry cache[kCacheN];
  thread_local int nfill = 0;

  for (int i = 0; i < nfill; ++i) {
    if (cache[i].kappa == kappa && cache[i].psi == psi) {
      const double z = cache[i].logZ;
      if (i > 0) {
        const Entry hit = cache[i];
        for (int k = i; k > 0; --k) {
          cache[k] = cache[k - 1];
        }
        cache[0] = hit;
      }
      return z;
    }
  }

  const double z = logNormalizer(kappa, psi);
  const int last = (nfill < kCacheN) ? nfill : (kCacheN - 1);
  for (int k = last; k > 0; --k) {
    cache[k] = cache[k - 1];
  }
  cache[0].kappa = kappa;
  cache[0].psi = psi;
  cache[0].logZ = z;
  if (nfill < kCacheN) {
    ++nfill;
  }
  return z;
}

inline double jonesPewseyLogDensity(double x, double mu, double kappa,
                                    double psi) {
  const double c = std::cos(wrapAngle(x) - wrapAngle(mu));
  return logUnnormalized(kappa, psi, c) - cachedLogNormalizer(kappa, psi);
}

}  // namespace jonespewsey
}  // namespace jags

#endif  /* JONESPEWSEY_MATH_H_ */
