#include "placer23b.h"
#include "log.h"

namespace replace
{
  /////////////////////////////////////////////////
  // Placer23b

  void Placer23b::printInfo() const
  {
    size_t maxFanout = 0;
    size_t sumFanout = 0;
    for (const Net23b* net : nets_)
    {
      maxFanout = std::max(net->pins().size(), maxFanout);
      sumFanout += net->pins().size();
    }

    LOG_INFO("NumInstances: {}", insts_.size());
    LOG_INFO("NumNets: {}", netStor_.size());
    LOG_INFO("NumPins: {}", sumFanout);
    LOG_INFO("MaxFanout: {}", maxFanout);
    LOG_INFO("AvgFanout: {}", sumFanout / (float)netStor_.size());

    LOG_INFO("DieCoreAreaLxLy: ({}, {})", topDie_.coreLx(), topDie_.coreLy());
    LOG_INFO("DieCoreAreaUxUy: ({}, {})", topDie_.coreUx(), topDie_.coreUy());
    int64_t coreArea =
        static_cast<int64_t>(topDie_.coreUx() - topDie_.coreLx()) *
        static_cast<int64_t>(topDie_.coreUy() - topDie_.coreLy());
    LOG_INFO("DieCoreArea: {}", coreArea);
    LOG_INFO("TopDieMaxUtil: {}", topUtil_);
    LOG_INFO("BottomDieMaxUtil: {}", bottomUtil_);

    LOG_INFO("TopDieTechnology Begin");
    topTech_->printInfo();
    LOG_INFO("TopDieTechnology End");
    LOG_INFO("BottomDieTechnology Begin");
    bottomTech_->printInfo();
    LOG_INFO("BottomDieTechnology End");

    LOG_INFO("TerminateSize: ({}, {})", termSizeX_, termSizeY_);
    LOG_INFO("TerminateSpacing: {}", termSpace_);
    LOG_INFO("TerminateCose: {}", termCost_);
  }
}
