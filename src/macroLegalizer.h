#ifndef __REPLACE_MACRO_LEGALIZER__
#define __REPLACE_MACRO_LEGALIZER__

#include <vector>
#include <memory>

namespace replace
{
  class PlacerBase;
  class NesterovBase;
  class Instance;
  class Die;
  class GCell;

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
    int overlapArea(GCell* , GCell*);
    int getCellMacroOverlap();
    int getMacrosOverlap();

  private:
    void postLegalize(const std::vector<Instance*> insts, Die* die);

  private:
    std::shared_ptr<PlacerBase> pb_;
    std::shared_ptr<NesterovBase> nb_;
    MacroLegalizerVars lgVars_;
  };
}

#endif