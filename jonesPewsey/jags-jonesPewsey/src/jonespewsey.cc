#include <module/Module.h>

#include "distributions/DJonesPewsey.h"
#include "functions/DJonesPewseyFun.h"

namespace jags {
namespace jonespewsey {

class JONESPEWSEYModule : public Module {
 public:
  JONESPEWSEYModule();
  ~JONESPEWSEYModule();
};

JONESPEWSEYModule::JONESPEWSEYModule() : Module("jonespewsey") {
  insert(new DJonesPewsey("dJonesPewsey"));
  insert(new DJonesPewsey("djonespewsey"));
  insert(new DJonesPewseyFun("dJonesPewsey"));
  insert(new DJonesPewseyFun("djonespewsey"));
}

JONESPEWSEYModule::~JONESPEWSEYModule() {
  std::vector<Distribution *> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  }
  std::vector<Function *> const &fvec = functions();
  for (unsigned int i = 0; i < fvec.size(); ++i) {
    delete fvec[i];
  }
}

}  // namespace jonespewsey
}  // namespace jags

jags::jonespewsey::JONESPEWSEYModule _jonespewsey_module;
