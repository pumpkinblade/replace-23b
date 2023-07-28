#ifndef __REPLACE_TERMINAL_MODIFIER__
#define __REPLACE_TERMINAL_MODIFIER__

#include <memory>

namespace replace
{
  class PlacerBase;

  class TerminalModifier
  {
  public:
    TerminalModifier() = default;
    ~TerminalModifier() = default;

    void setPlacerBase(std::shared_ptr<PlacerBase> pb) { pb_ = pb; }

    void modify();
    void recover();

  private:
    std::shared_ptr<PlacerBase> pb_;
  };  
}

#endif