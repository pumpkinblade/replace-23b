#ifndef __REPLACE_MACRO_LEGALIZER__
#define __REPLACE_MACRO_LEGALIZER__

#include <vector>
#include <memory>

namespace replace
{
  class PlacerBase;
  class Instance;
  class Die;

  class MacroLegalizerVars
  {
  public:
    MacroLegalizerVars();

  public:
    int maxPostLegalizeIter;
  };

  class MacroLegalizer
  {
  public:
    MacroLegalizer(MacroLegalizerVars lgVars, std::shared_ptr<PlacerBase> pb);

    void doLegalization();

  private:
    void postLegalize(const std::vector<Instance*> insts, Die* die);

  private:
    std::shared_ptr<PlacerBase> pb_;
    MacroLegalizerVars lgVars_;
  };
}

#endif