#ifndef DJONESPEWSEY_FUN_H_
#define DJONESPEWSEY_FUN_H_

#include <string>

#include <function/ScalarFunction.h>

namespace jags {
namespace jonespewsey {

/**
 * Density function matching R's djonespewsey(x, mu, kappa, psi).
 *
 * Usage: y <- dJonesPewsey(x, mu, kappa, psi)
 */
class DJonesPewseyFun : public ScalarFunction {
 public:
  explicit DJonesPewseyFun(std::string const &name = "dJonesPewsey");
  double evaluate(std::vector<double const *> const &args) const;
  bool checkParameterValue(std::vector<double const *> const &args) const;
};

}  // namespace jonespewsey
}  // namespace jags

#endif  /* DJONESPEWSEY_FUN_H_ */
