#include "DJonesPewseyFun.h"

#include <cmath>

#include "../JonesPewseyMath.h"

namespace jags {
namespace jonespewsey {

DJonesPewseyFun::DJonesPewseyFun(std::string const &name)
    : ScalarFunction(name, 4) {}

double DJonesPewseyFun::evaluate(
    std::vector<double const *> const &args) const {
  const double ld = jonesPewseyLogDensity(*args[0], *args[1], *args[2],
                                          *args[3]);
  if (!std::isfinite(ld)) {
    return 0.0;
  }
  return std::exp(ld);
}

bool DJonesPewseyFun::checkParameterValue(
    std::vector<double const *> const &args) const {
  return std::isfinite(*args[0]) && std::isfinite(*args[1]) &&
         std::isfinite(*args[2]) && std::isfinite(*args[3]) &&
         *args[2] >= 0.0;
}

}  // namespace jonespewsey
}  // namespace jags
