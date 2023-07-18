#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "placerBase.h"
#include "parser.h"
#include "log.h"
#include "replace.h"
#include "plot.h"

using namespace replace;
using namespace TCLAP;
using namespace std;

int main(int argc, const char *argv[])
{
  Log::Init();

  PlotVars vars;
  vars.minLength = 1000;
  vars.xMargin = 30;
  vars.yMargin = 30;
  Plot::init(vars);

  string lefFilename;
  string defFilename;
  string mode;
  string txtFilename;
  float  targetDensity;

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "0.1(alpha)");

    // Define a value argument and add it to the command line.
    ValueArg<string> lefArg("l", "lef", "path to lef file", false, "none", "string");
    ValueArg<string> defArg("d", "def", "path to def file", false, "none", "string");
    ValueArg<string> modeArg("m", "mode", "lefdef/23b", false, "lefdef", "string");
    ValueArg<string> txtArg("b", "txt23b", "path to 23b text file", false, "none", "string");
    ValueArg<float> densityArg("D", "density", "target density", false, 1.0, "float");

    cmd.add(lefArg);
    cmd.add(defArg);
    cmd.add(modeArg);
    cmd.add(txtArg);
    cmd.add(densityArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    lefFilename = lefArg.getValue();
    defFilename = defArg.getValue();
    mode = modeArg.getValue();
    txtFilename = txtArg.getValue();
    targetDensity = densityArg.getValue();
  }
  catch (ArgException &e) // catch any exceptions
  {
    LOG_ERROR("TCLAP Error: {} for arg {}", e.error(), e.argId());
    return 1;
  }

  if (mode == "lefdef")
  {
    LOG_TRACE("Parse Lef/Def Begin");
    std::shared_ptr<PlacerBase> pb = Parser::lefdefToPlacerBase(lefFilename, defFilename);
    LOG_TRACE("Parse Lef/Def End");
    pb->printInfo();

    Replace rp(targetDensity);
    rp.setPlacerBase(pb);
    rp.doInitialPlace();
    rp.doNesterovPlace();
    rp.doAbacusLegalization();
  }
  else if (mode == "23b")
  {
    LOG_TRACE("Parse 23b Text File Begin");
    std::shared_ptr<PlacerBase> pb = Parser::txtToPlacerBase(txtFilename);
    pb->printInfo();
    LOG_TRACE("Parse 23b Text File End");

    Replace rp(targetDensity);
    rp.setPlacerBase(pb);
    rp.doInitialPlace();
    rp.doNesterovPlace();
  }

  return 0;
}
