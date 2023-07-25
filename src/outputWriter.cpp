#include "outputWriter.h"
#include "placerBase.h"
#include "log.h"
#include <fstream>

namespace replace
{
  void writeInstance(std::ofstream& out, const Instance* inst)
  {
    out << "Inst " << inst->name() << " " 
        << inst->lx() << " " << inst->ly() << " R0\n";
  }

  void writeTerminal(std::ofstream& out, const Instance* term)
  {
    out << "Terminal " << term->name() << " " 
        << term->cx() << " " << term->cy() << "\n";
  }

  void OutputWriter::write(PlacerBase* pb, std::string& filename)
  {
    std::ofstream out(filename);
    if(out.bad())
    {
      LOG_ERROR("Couldn't open file `{}`", filename);
      return;
    }

    // top die
    const Die* topDie = pb->die("top");
    if(topDie == nullptr)
    {
      LOG_ERROR("Couldn't find top die");
      return;
    }
    out << "TopDiePlacement " << topDie->insts().size() << '\n';
    for(const Instance* inst : topDie->insts())
    {
      writeInstance(out, inst);
    }
    out << std::endl;

    // bottom die
    const Die* bottomDie = pb->die("bottom");
    if(bottomDie == nullptr)
    {
      LOG_ERROR("Couldn't find bottom die");
      return;
    }
    out << "BottomDiePlacement " << bottomDie->insts().size() << '\n';
    for(const Instance* inst : bottomDie->insts())
    {
      writeInstance(out, inst);
    }
    out << std::endl;

    // terminal die
    const Die* terminalDie = pb->die("terminal");
    if(terminalDie == nullptr)
    {
      LOG_ERROR("Couldn't find top die");
      return;
    }
    out << "NumTerminals " << terminalDie->insts().size() << '\n';
    for(const Instance* term : terminalDie->insts())
    {
      writeTerminal(out, term);
    }
    out << std::endl;
  }
}