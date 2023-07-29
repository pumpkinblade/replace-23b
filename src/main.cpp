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

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "beta");

    // Define a value argument and add it to the command line.
    ValueArg<string> inputArg("i", "input", "path to input file", false, "none", "string");
    ValueArg<string> outputArg("o", "output", "path to output file", false, "none", "string");
    cmd.add(inputArg);
    cmd.add(outputArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    inputFilename = inputArg.getValue();
    outputFilename = outputArg.getValue();
  }
  catch (ArgException &e) // catch any exceptions
  {
    LOG_CRITICAL("TCLAP Error: {} for arg {}", e.error(), e.argId());
    return 1;
  }

  LOG_TRACE("Parse 23b Text File Begin");
  std::shared_ptr<PlacerBase> pb = Parser::txtToPlacerBase(inputFilename);
  pb->printDebugInfo();
  LOG_TRACE("Parse 23b Text File End");

  // then we do partition
  Partitioner partitioner(1.0f);
  partitioner.partitioning2(pb);
  Plot::plot(pb.get(), "./plot/cell", "part");
  
  // then we do optimization
  TerminalModifier tm;
  tm.setPlacerBase(pb);
  tm.modify();
  Replace rp(1.0f);
  rp.setPlacerBase(pb);
  rp.doInitialPlace();
  LOG_INFO("HPWL after ip: {}", pb->hpwl());
  Plot::plot(pb.get(), "./plot/cell", "ip");
  rp.doNesterovPlace();
  LOG_INFO("HPWL after gp: {}", pb->hpwl());
  Plot::plot(pb.get(), "./plot/cell", "gp");
  rp.doMacroLegalization();
  LOG_INFO("HPWL after mlg: {}", pb->hpwl());
  Plot::plot(pb.get(), "./plot/cell", "mlg");
  rp.doAbacusLegalization();
  LOG_INFO("HPWL after clg: {}", pb->hpwl());
  Plot::plot(pb.get(), "./plot/cell", "clg");
  tm.recover();

  LOG_INFO("Result HPWL: {}", pb->hpwl());
  Plot::plot(pb.get(), "./plot/cell", "res");

  OutputWriter::write(pb.get(), outputFilename);

  return 0;
}
 
