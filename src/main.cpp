#include <iostream>
#include <string>
#include <tclap/CmdLine.h>
#include <omp.h>

#include "placerBase.h"
#include "parser.h"
#include "log.h"
#include "replace.h"
#include "plot.h"
#include "partitioner.h"
#include "terminalModifier.h"
#include "outputWriter.h"

using namespace replace;
using namespace TCLAP;
using namespace std;

int main(int argc, const char *argv[])
{
  Log::Init();
#if defined(DEBUG) || defined(_DEBUG)
  Log::getLogger()->set_level(spdlog::level::trace);
#else
  Log::getLogger()->set_level(spdlog::level::info);
#endif

  PlotVars vars;
  vars.minLength = 1000;
  vars.xMargin = 30;
  vars.yMargin = 30;
  Plot::init(vars);

  omp_set_num_threads(16);

  string inputFilename;
  string outputFilename;
  string mode;

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "beta");

    // Define a value argument and add it to the command line.
    ValueArg<string> inputArg("i", "input", "path to input file", false, "none", "string");
    ValueArg<string> outputArg("o", "output", "path to output file", false, "none", "string");
    ValueArg<string> modeArg("m", "mode", "lefdef/23b", false, "lefdef", "string");
    cmd.add(inputArg);
    cmd.add(outputArg);
    cmd.add(modeArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    inputFilename = inputArg.getValue();
    outputFilename = outputArg.getValue();
    mode = modeArg.getValue();
  }
  catch (ArgException &e) // catch any exceptions
  {
    LOG_CRITICAL("TCLAP Error: {} for arg {}", e.error(), e.argId());
    return 1;
  }

  // if (mode == "lefdef")
  // {
  //   LOG_TRACE("Parse Lef/Def Begin");
  //   std::shared_ptr<PlacerBase> pb = Parser::lefdefToPlacerBase(lefFilename, defFilename);
  //   LOG_TRACE("Parse Lef/Def End");
  //   pb->printDebugInfo();

  //   Replace rp(targetDensity);
  //   rp.setPlacerBase(pb);
  //   rp.doInitialPlace();
  //   rp.doNesterovPlace();
  //   rp.doAbacusLegalization();
  // }
  // else if (mode == "23b")
  if(mode == "23b")
  {
    LOG_TRACE("Parse 23b Text File Begin");
    std::shared_ptr<PlacerBase> pb = Parser::txtToPlacerBase(inputFilename);
    pb->printDebugInfo();
    LOG_TRACE("Parse 23b Text File End");

    // then we do partition
    Partitioner partitioner(1.0);
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

    rp.doNesterovPlace("finalgp");
    Plot::plot(pb.get(), "./plot/cell", "after_finalgp");
    rp.doAbacusLegalization();
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
    std::shared_ptr<PlacerBase> pb = Parser::txtToPlacerBase(inputFilename);
    pb->printDebugInfo();
    LOG_TRACE("Parse 23b Text File End");

    // then we do partition
    Partitioner partitioner(1.0);
    partitioner.hmetistest(pb);

    // Replace rp(targetDensity);
    // rp.setPlacerBase(pb);
    // rp.doInitialPlace();
    // rp.doNesterovPlace("pregp");
    // rp.setTargetDensity(1.0);
    // rp.doNesterovPlace("postgp");
  }


  return 0;
}
 
