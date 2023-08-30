#include <iostream>
#include <string>
#include <unordered_set>
#include <tclap/CmdLine.h>
#include <spdlog/sinks/basic_file_sink.h>

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
  string outputFilename;
  bool useTheta;
  bool useLocalDensity;
  string statFilename;

  // Wrap everything in a try block.  Do this every time, because exceptions will be thrown for problems.
  try
  {
    CmdLine cmd("Command description message", ' ', "beta");

    // Define a value argument and add it to the command line.
    ValueArg<string> txtArg("i", "input", "path to input file", false, "none", "string");
    ValueArg<string> outputArg("o", "output", "path to output file", false, "output.txt", "string");
    SwitchArg thetaArg("R", "rotate", "enable rotation", false);
    SwitchArg localArg("L", "local", "enable local density", false);
    ValueArg<string> statArg("s", "statistic", "path to statistic file", false, "statistic.log", "string");

    cmd.add(txtArg);
    cmd.add(outputArg);
    cmd.add(thetaArg);
    cmd.add(localArg);
    cmd.add(statArg);

    // Parse the args.
    cmd.parse(argc, argv);

    // Get the value parsed by each arg.
    txtFilename = txtArg.getValue();
    outputFilename = outputArg.getValue();
    useTheta = thetaArg.getValue();
    useLocalDensity = localArg.getValue();
    statFilename = statArg.getValue();
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
  Partitioner partitioner();
  // double x=partitioner.getAverageTechRatio(pb);
  // LOG_INFO("Average Tech Ratio: {}", x);
  // partitioner.partitioning2(pb);
  // partitioner.partitionInstance(pb);
  // partitioner.mtPartitionInstance(pb);
  Partitioner partitioner = Partitioner();
  // partitioner.partitioning2(pb);
#if defined(WIN32) || defined(_WIN32)
  partitioner.partitionInstance(pb);
#else
  partitioner.mtPartitionInstance2(pb);
#endif
  Plot::plot(pb.get(), "./plot/cell", "after_partition");

  // then we do optimization
  Replace rp(1.0);
  rp.setPlacerBase(pb);
  rp.modifyTerminal();
  rp.setInitialPlaceMinIter(5);
  rp.doInitialPlace();
  rp.setNesterovUseTheta(useTheta);
  rp.setNesterovPlaceUseLocalDensity(useLocalDensity);
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
  rp.setNesterovUseTheta(false);
  rp.setNesterovPlaceUseLocalDensity(useLocalDensity);
  rp.doNesterovPlace("finalgp");
  Plot::plot(pb.get(), "./plot/cell", "after_finalgp");
  rp.doAbacusLegalization();
  Plot::plot(pb.get(), "./plot/cell", "after_clg");
  rp.recoverTerminal();

  int64_t hpwl = pb->hpwl();
  LOG_INFO("Result HPWL: {}", hpwl);
  Plot::plot(pb.get(), "./plot/cell", "result");

  OutputWriter::write(pb.get(), outputFilename);

  {
    // find the largest net
    auto& nets = pb->nets();
    std::sort(nets.begin(), nets.end(), [](const Net *net1, const Net *net2)
                                        { return net1->hpwl() > net2->hpwl(); });
    for(int i = 0; i < 10; i++)
    {
      const Net* lnet = nets[i];
      // check if the net cross dies
      const Instance* term = nullptr;
      for(const Pin* pin : lnet->pins())
      {
        if(pin->instance()->name().front() == 'N')
        {
          term = pin->instance();
          break;
        }
      }
      if(term)
      {
        Plot::plotCNetR1(pb.get(), term, "./plot/net", "net" + std::to_string(i));
      }
      else
      {
        const Instance* inst = lnet->pins().front()->instance();
        if(std::find(pb->die("top")->insts().begin(), pb->die("top")->insts().end(), inst) != pb->die("top")->insts().end())
        {
          Plot::plotNetR1(pb.get(), pb->die("top"), lnet, "./plot/net", "net" + std::to_string(i));
        }
        else
        {
          Plot::plotNetR1(pb.get(), pb->die("bottom"), lnet, "./plot/net", "net" + std::to_string(i));
        }
      }
    }

    auto statSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(statFilename, true);
    statSink->set_pattern("%v");
    auto statLogger = std::make_shared<spdlog::logger>("stat", statSink);
    spdlog::register_logger(statLogger);

    // net hpwl
    double maxNetHpwl, minNetHpwl, avgNetHpwl, totalHpwl;
    double maxTopNetHpwl, minTopNetHpwl, avgTopNetHpwl;
    double maxBotNetHpwl, minBotNetHpwl, avgBotNetHpwl;
    maxNetHpwl = maxTopNetHpwl = maxBotNetHpwl = std::numeric_limits<double>::lowest();
    minNetHpwl = minTopNetHpwl = minBotNetHpwl = std::numeric_limits<double>::max();
    totalHpwl = avgNetHpwl = avgTopNetHpwl = avgBotNetHpwl = 0;
    for (const Net *net : pb->nets())
    {
      double hpwl = static_cast<double>(net->hpwl());
      maxNetHpwl = std::max(hpwl, maxNetHpwl);
      minNetHpwl = std::min(hpwl, minNetHpwl);
      totalHpwl += hpwl;
    }
    avgNetHpwl = totalHpwl / static_cast<double>(pb->nets().size());
    for (const Instance *inst : pb->die("terminal")->insts())
    {
      const Pin *topPin = inst->pins()[0];
      const Pin *botPin = inst->pins()[1];
      const Net *topNet = topPin->net();
      const Net *botNet = botPin->net();
      double topHpwl = static_cast<double>(topNet->hpwl());
      double botHpwl = static_cast<double>(botNet->hpwl());
      maxTopNetHpwl = std::max(topHpwl, maxTopNetHpwl);
      minTopNetHpwl = std::min(topHpwl, minTopNetHpwl);
      avgTopNetHpwl += topHpwl;
      maxBotNetHpwl = std::max(botHpwl, maxBotNetHpwl);
      minBotNetHpwl = std::min(botHpwl, minBotNetHpwl);
      avgBotNetHpwl += botHpwl;
    }
    avgTopNetHpwl /= static_cast<double>(pb->die("terminal")->insts().size());
    avgBotNetHpwl /= static_cast<double>(pb->die("terminal")->insts().size());

    statLogger->info("case {}", txtFilename);
    statLogger->info("AvgNetHpwl {}", avgNetHpwl);
    statLogger->info("MaxNetHpwl {}", maxNetHpwl);
    statLogger->info("MinNetHpwl {}", minNetHpwl);
    statLogger->info("AvgTopNetHpwl {}", avgTopNetHpwl);
    statLogger->info("MaxTopNetHpwl {}", maxTopNetHpwl);
    statLogger->info("MinTopNetHpwl {}", minTopNetHpwl);
    statLogger->info("AvgBotNetHpwl {}", avgBotNetHpwl);
    statLogger->info("MaxBotNetHpwl {}", maxBotNetHpwl);
    statLogger->info("MinBotNetHpwl {}", minBotNetHpwl);
    statLogger->info("TotalHpwl {}", totalHpwl);
    statLogger->info("TerminalCost {}", pb->terminalCost());
    statLogger->info("TerminalCount {}", pb->die("terminal")->insts().size());

    statLogger->info("================= Net =======================");
    statLogger->info("hpwl, num_insts, num_macros");
    for(const Net* net : pb->nets())
    {
      std::unordered_set<const Instance*> netInstances;
      std::unordered_set<const Instance*> netMacros;
      for(const Pin* pin : net->pins())
      {
        const Instance* inst = pin->instance();
        netInstances.insert(inst);
        if(inst->isMacro())
          netMacros.insert(inst);
      }
      statLogger->info("{}, {}, {}", net->hpwl(), netInstances.size(), netMacros.size());
    }

    statLogger->info("================= Net With Terminal =======================");
    statLogger->info("hpwl_top, hpwl_bot, num_insts_top, num_insts_bot, num_macros_top, num_macros_bot");
    for(const Instance* term : pb->die("terminal")->insts())
    {
      const Net* net0 = term->pins()[0]->net();
      const Net* net1 = term->pins()[1]->net();
      std::unordered_set<const Instance*> topNetInstances;
      std::unordered_set<const Instance*> topNetMacros;
      std::unordered_set<const Instance*> botNetInstances;
      std::unordered_set<const Instance*> botNetMacros;
      for(const Pin* pin : net0->pins())
      {
        const Instance* inst = pin->instance();
        topNetInstances.insert(inst);
        if (inst->isMacro())
          topNetMacros.insert(inst);
      }
      for(const Pin* pin : net1->pins())
      {
        const Instance* inst = pin->instance();
        botNetInstances.insert(inst);
        if (inst->isMacro())
          botNetMacros.insert(inst);
      }
      statLogger->info("{}, {}, {}, {}, {}, {}", net0->hpwl(), net1->hpwl(),
                       topNetInstances.size(), botNetInstances.size(),
                       topNetMacros.size(), botNetMacros.size());
    }
  }

  return 0;
}
