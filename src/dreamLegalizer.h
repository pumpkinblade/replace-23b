#ifndef __REPLACE_DREAM_LEGALIZER__
#define __REPLACE_DREAM_LEGALIZER__

#include <memory>

namespace replace
{
  class PlacerBase;

  class DreamLegalizerVars
  {
  public:
    enum { One, Area, NumPins } weightOpt;

    DreamLegalizerVars();
  };

  class DreamLegalizer
  {
  public:
    DreamLegalizer() = default;
    DreamLegalizer(DreamLegalizerVars vars, std::shared_ptr<PlacerBase> pb);
    ~DreamLegalizer() = default;

    void doLegalization();
  
  private:
    std::shared_ptr<PlacerBase> pb_;
    DreamLegalizerVars vars_;
  };
}

#endif