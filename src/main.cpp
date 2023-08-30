#include <iostream>
#include <string>
#include <tclap/CmdLine.h>

#include "placerBase.h"
#include "parser.h"
#include "log.h"
#include "replace.h"
#include "plot.h"
#include "partitioner.h"
#include "outputWriter.h"
#include "abaxLegalizer.h"

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

  string txtFilename;
  float targetDensity;
  string outputFilename;

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "beta");

    // Define a value argument and add it to the command line.
    ValueArg<string> txtArg("i", "input", "path to input file", false, "none", "string");
    ValueArg<float> densityArg("D", "density", "target density", false, 1.0, "float");
    ValueArg<string> outputArg("o", "output", "path to output file", false, "output.txt", "string");

    cmd.add(txtArg);
    cmd.add(densityArg);
    cmd.add(outputArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    txtFilename = txtArg.getValue();
    targetDensity = densityArg.getValue();
    outputFilename = outputArg.getValue();
  }
  catch (ArgException &e) // catch any exceptions
  {
    LOG_ERROR("TCLAP Error: {} for arg {}", e.error(), e.argId());
    return 1;
  }

  LOG_TRACE("Parse 23b Text File Begin");
  std::shared_ptr<PlacerBase> pb = Parser::txtToPlacerBase(txtFilename);
  pb->printDebugInfo();
  LOG_TRACE("Parse 23b Text File End");

  // then we do partition
  Partitioner partitioner(targetDensity);
  // partitioner.partitioning2(pb);
  partitioner.partitionInstance(pb);
  // partitioner.mtPartitionInstance(pb);
  Plot::plot(pb.get(), "./plot/cell", "after_partition");

  // then we do optimization
  Replace rp(1.0);
  rp.setPlacerBase(pb);
  rp.modifyTerminal();
  rp.setInitialPlaceMinIter(5);
  rp.doInitialPlace();
  rp.setNesterovUseTheta(true);
  rp.setNesterovPlaceUseLocalDensity(true);
  rp.doNesterovPlace("postgp");
  Plot::plot(pb.get(), "./plot/cell", "after_postgp");

  rp.setMacroPostMaxIter(100000);
  rp.doMacroLegalization();
  Plot::plot(pb.get(), "./plot/cell", "after_mlg");

  // fix macros
  for (Instance* inst : pb->insts())
  {
    if (inst->isMacro())
      inst->setFixed(true);
  }

  rp.doInitialPlace();
  rp.setNesterovPlaceUseLocalDensity(true);
  rp.doNesterovPlace("finalgp");
  Plot::plot(pb.get(), "./plot/cell", "after_finalgp");
  rp.doAbacusLegalization();
  Plot::plot(pb.get(), "./plot/cell", "after_clg");
  rp.recoverTerminal();

  int64_t hpwl = pb->hpwl();
  LOG_INFO("Result HPWL: {}", hpwl);
  Plot::plot(pb.get(), "./plot/cell", "result");

  OutputWriter::write(pb.get(), outputFilename);

  return 0;
}
