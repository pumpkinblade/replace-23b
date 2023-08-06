#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "placerBase.h"
#include "parser.h"
#include "log.h"
#include "replace.h"
#include "plot.h"
#include "partitioner.h"
#include "terminalModifier.h"
#include "outputWriter.h"
#include "dreamLegalizer.h"

using namespace replace;
using namespace TCLAP;
using namespace std;

int main(int argc, const char *argv[])
{
  Log::Init();

  PlotVars vars;
  vars.minLength = 2000;
  vars.xMargin = 30;
  vars.yMargin = 30;
  Plot::init(vars);

  string lefFilename;
  string defFilename;
  string mode;
  string txtFilename;
  float  targetDensity;
  string outputFilename;

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "beta");

    // Define a value argument and add it to the command line.
    ValueArg<string> lefArg("l", "lef", "path to lef file", false, "none", "string");
    ValueArg<string> defArg("d", "def", "path to def file", false, "none", "string");
    ValueArg<string> modeArg("m", "mode", "lefdef/23b", false, "lefdef", "string");
    ValueArg<string> txtArg("b", "txt23b", "path to 23b text file", false, "none", "string");
    ValueArg<float> densityArg("D", "density", "target density", false, 1.0, "float");
    ValueArg<string> outputArg("o", "output", "path to output file", false, "output.txt", "string");

    cmd.add(lefArg);
    cmd.add(defArg);
    cmd.add(modeArg);
    cmd.add(txtArg);
    cmd.add(densityArg);
    cmd.add(outputArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    lefFilename = lefArg.getValue();
    defFilename = defArg.getValue();
    mode = modeArg.getValue();
    txtFilename = txtArg.getValue();
    targetDensity = densityArg.getValue();
    outputFilename = outputArg.getValue();
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
    pb->printDebugInfo();

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
    pb->printDebugInfo();
    LOG_TRACE("Parse 23b Text File End");

    // then we do partition
    Partitioner partitioner(targetDensity);
    partitioner.partitioning2(pb);
    Plot::plot(pb.get(), "./plot/cell", "after_partition");
    
    // then we do optimization
    TerminalModifier tm;
    tm.setPlacerBase(pb);
    tm.modify();
    Replace rp(1.0);
    rp.setPlacerBase(pb);
    rp.doInitialPlace();
    rp.doNesterovPlace("postgp");
    Plot::plot(pb.get(), "./plot/cell", "after_postgp");
    rp.doMacroLegalization();
    Plot::plot(pb.get(), "./plot/cell", "after_mlg");

    // fix macros
    {
      for(Instance* inst : pb->insts())
      {
        if(inst->isMacro())
          inst->setFixed(true);
      }
    }
    
    rp.setTargetOverflow(0.05f);
    rp.setIncrementalPlaceMode(true);
    rp.doInitialPlace();
    rp.doNesterovPlace("finalgp");
    Plot::plot(pb.get(), "./plot/cell", "after_finalgp");
    rp.doAbacusLegalization();
    //{
    //  DreamLegalizer dlg(DreamLegalizerVars(), pb);
    //  dlg.doLegalization();
    //}
    Plot::plot(pb.get(), "./plot/cell", "after_clg");
    tm.recover();

    int64_t hpwl = pb->hpwl();
    LOG_INFO("Result HPWL: {}", hpwl);
    Plot::plot(pb.get(), "./plot/cell", "result");

    OutputWriter::write(pb.get(), outputFilename);
  }
  else if (mode == "latest")
  {
    LOG_TRACE("Parse 23b Text File Begin");
    std::shared_ptr<PlacerBase> pb = Parser::txtToPlacerBase(txtFilename);
    pb->printDebugInfo();
    LOG_TRACE("Parse 23b Text File End");

    Replace rp(targetDensity);
    rp.setPlacerBase(pb);
    rp.doInitialPlace();
    rp.doNesterovPlace("pregp");
    rp.setTargetDensity(1.0);
    rp.doNesterovPlace("postgp");
  }


  return 0;
}
