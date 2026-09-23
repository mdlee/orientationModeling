#ifndef DJONESPEWSEY_H_
#define DJONESPEWSEY_H_

#include <string>

#include <distribution/ScalarDist.h>

namespace jags {
namespace jonespewsey {

/**
 * Jones-Pewsey symmetric circular distribution on [0, 2*pi).
 *
 * Parameters:
 *   mu    - mean direction; any real value, wrapped internally
 *   kappa - concentration, kappa >= 0
 *   psi   - real shape (continuous). psi -> 0 is von Mises,
 *           psi = -1 wrapped Cauchy, psi = 1 cardioid.
 *
 * Density:
 *   f(x) = [cosh(kappa*psi) + sinh(kappa*psi)*cos(x-mu)]^(1/psi)
 *          / (2*pi * P_{1/psi}(cosh(kappa*psi)))
 *
 * Reference:
 *   Jones, M. C. & Pewsey, A. (2005). A family of symmetric
 *   distributions on the circle. JASA 100, 1422-1428.
 */
class DJonesPewsey : public ScalarDist {
 public:
  explicit DJonesPewsey(std::string const &name = "dJonesPewsey");

  double logDensity(double x, PDFType type,
                    std::vector<double const *> const &parameters,
                    double const *lower, double const *upper) const;
  double randomSample(std::vector<double const *> const &parameters,
                      double const *lower, double const *upper,
                      RNG *rng) const;
  double typicalValue(std::vector<double const *> const &parameters,
                      double const *lower, double const *upper) const;
  bool checkParameterValue(std::vector<double const *> const &parameters) const;

  double l(std::vector<double const *> const &parameters) const;
  double u(std::vector<double const *> const &parameters) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

}  // namespace jonespewsey
}  // namespace jags

#endif  /* DJONESPEWSEY_H_ */
