#include "terminalModifier.h"
#include "placerBase.h"

namespace replace
{
  void TerminalModifier::modify()
  {
    // terminal resizing with spacing
    int spaceLeft = pb_->terminalSpacing();
    int spaceRight = pb_->terminalSpacing() - spaceLeft;
    for(Instance* inst : pb_->die("terminal")->insts())
    {
      inst->setBox(
        inst->lx() - spaceLeft,
        inst->ly() - spaceLeft,
        inst->ux() + spaceRight,
        inst->uy() + spaceRight
      );
    }

    // terminal die core area
    int coreLx = pb_->die("terminal")->coreLx() + spaceRight;
    int coreLy = pb_->die("terminal")->coreLy() + spaceRight;
    int coreUx = pb_->die("terminal")->coreUx() - spaceLeft;
    int coreUy = pb_->die("terminal")->coreUy() - spaceLeft;
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
    int spaceLeft = pb_->terminalSpacing();
    int spaceRight = pb_->terminalSpacing() - spaceLeft;
    for(Instance* inst : pb_->die("terminal")->insts())
    {
      inst->setBox(
        inst->lx() + spaceLeft,
        inst->ly() + spaceLeft,
        inst->ux() - spaceRight,
        inst->uy() - spaceRight
      );
    }

    // terminal die core area
    int coreLx = pb_->die("terminal")->coreLx() - spaceRight;
    int coreLy = pb_->die("terminal")->coreLy() - spaceRight;
    int coreUx = pb_->die("terminal")->coreUx() + spaceLeft;
    int coreUy = pb_->die("terminal")->coreUy() + spaceLeft;
    pb_->die("terminal")->setCoreBox(coreLx, coreLy, coreUx, coreUy);
  }
}