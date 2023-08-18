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
#include "abaxLegalizer.h"

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

  string txtFilename;
  float targetDensity;
  string outputFilename;
  string mode;

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "beta");

    // Define a value argument and add it to the command line.
    ValueArg<string> txtArg("i", "input", "path to input file", false, "none", "string");
    ValueArg<float> densityArg("D", "density", "target density", false, 1.0, "float");
    ValueArg<string> outputArg("o", "output", "path to output file", false, "output.txt", "string");
    ValueArg<string> modeArg("m", "mode", "mode", false, "23b/latest", "string");

    cmd.add(txtArg);
    cmd.add(densityArg);
    cmd.add(outputArg);
    cmd.add(modeArg);


    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    txtFilename = txtArg.getValue();
    targetDensity = densityArg.getValue();
    outputFilename = outputArg.getValue();
    mode = modeArg.getValue();
    if (mode == "23b")
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
      //rp.doAbacusLegalization();
      //{
      //  DreamLegalizer dlg(DreamLegalizerVars(), pb);
      //  dlg.doLegalization();
      //}
      {
        AbaxLegalizer abax(pb);
        abax.doLegalization();
      }
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

      // then we do partition
      Partitioner partitioner(1.0);
      // partitioner.hmetistest(pb);
      // partitioner.do_run_kahypar();
      partitioner.partitionInstance(pb);
      // Plot::plot(pb.get(), "./plot/cell", "after_partition");
      // // then we do optimization
      // LOG_INFO("optimization");
      // TerminalModifier tm;
      // tm.setPlacerBase(pb);
      // tm.modify();
      // Replace rp(1.0);
      // rp.setPlacerBase(pb);
      // rp.doInitialPlace();
      // rp.doNesterovPlace("postgp");
      // Plot::plot(pb.get(), "./plot/cell", "after_postgp");
      // rp.doMacroLegalization();
      // Plot::plot(pb.get(), "./plot/cell", "after_mlg");

      // // fix macros
      // {
      //   for(Instance* inst : pb->insts())
      //   {
      //     if(inst->isMacro())
      //       inst->setFixed(true);
      //   }
      // }

      // LOG_INFO("finalgp");
      // rp.doNesterovPlace("finalgp");
      // Plot::plot(pb.get(), "./plot/cell", "after_finalgp");

      // // rp.doAbacusLegalization();
      // {
      //   AbaxLegalizer abax(pb);
      //   abax.doLegalization();
      // }
      // Plot::plot(pb.get(), "./plot/cell", "after_clg");
      // tm.recover();
      // int64_t hpwl = pb->hpwl();
      // LOG_INFO("Result HPWL: {}", hpwl);
      // Plot::plot(pb.get(), "./plot/cell", "result");

      // OutputWriter::write(pb.get(), outputFilename);
    }
    else
    {
      LOG_ERROR("Unknown mode: {}", mode);
    }
  }
  catch (ArgException &e) // catch any exceptions
  {
    LOG_ERROR("TCLAP Error: {} for arg {}", e.error(), e.argId());
    return 1;
  }


  return 0;
}
