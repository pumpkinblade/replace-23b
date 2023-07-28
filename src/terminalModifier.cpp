#include "terminalModifier.h"
#include "placerBase.h"

namespace replace
{
  void TerminalModifier::modify()
  {
    // terminal resizing with spacing
    int modSizeX = pb_->terminalSizeX() + pb_->terminalSpacing();
    int modSizeY = pb_->terminalSizeX() + pb_->terminalSpacing();
    for(Instance* inst : pb_->die("terminal")->insts())
    {
      int cx = inst->cx();
      int cy = inst->cy();
      inst->setSize(modSizeX, modSizeY);
      inst->setCenterLocation(cx, cy);
    }

    // terminal die core area
    int coreLx = pb_->die("terminal")->coreLx() + pb_->terminalSpacing() / 2;
    int coreLy = pb_->die("terminal")->coreLy() + pb_->terminalSpacing() / 2;
    int coreUx = pb_->die("terminal")->coreUx() - pb_->terminalSpacing() / 2;
    int coreUy = pb_->die("terminal")->coreUy() - pb_->terminalSpacing() / 2;
    pb_->die("terminal")->setCoreBox(coreLx, coreLy, coreUx, coreUy);

    // set row Params
    int rowWidth = pb_->die("terminal")->coreDx();
    int rowHeight = pb_->terminalSizeX() + pb_->terminalSpacing();
    int repeatCount = pb_->die("terminal")->coreDy() / rowHeight;
    int rowStartX = pb_->die("terminal")->coreLx();
    int rowStartY = pb_->die("terminal")->coreLy();
    pb_->die("terminal")->setRowParams(rowStartX, rowStartY, rowWidth, rowHeight, repeatCount);
  }

  void TerminalModifier::recover()
  {
    int origSizeX = pb_->terminalSizeX();
    int origSizeY = pb_->terminalSizeY();
    for(Instance* inst : pb_->die("terminal")->insts())
    {
      int cx = inst->cx();
      int cy = inst->cy();
      inst->setSize(origSizeX, origSizeY);
      inst->setCenterLocation(cx, cy);
    }

    // terminal die core area
    int coreLx = pb_->die("terminal")->coreLx() - pb_->terminalSpacing() / 2;
    int coreLy = pb_->die("terminal")->coreLy() - pb_->terminalSpacing() / 2;
    int coreUx = pb_->die("terminal")->coreUx() + pb_->terminalSpacing() / 2;
    int coreUy = pb_->die("terminal")->coreUy() + pb_->terminalSpacing() / 2;
    pb_->die("terminal")->setCoreBox(coreLx, coreLy, coreUx, coreUy);
  }
}