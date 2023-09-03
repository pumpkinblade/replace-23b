#ifndef __PARTITIONER__
#define __PARTITIONER__

#include <vector>
#include "technology.h"
#include "placerBase.h"

namespace replace
{
  class Replace;

  class Partitioner
  {
  public:
    Partitioner() = default;
    ~Partitioner() = default;

    void partitioning2(std::shared_ptr<PlacerBase> pb_);

    void partitionInstance(std::shared_ptr<PlacerBase> &pb_);
    void partitionInstance1(std::shared_ptr<PlacerBase> &pb_);
#if !defined(WIN32) && !defined(_WIN32)
    void mtPartitionInstance(std::shared_ptr<PlacerBase> &pb_);

    void mtPartitionInstance2(std::shared_ptr<PlacerBase> &pb_);
    void partitionTest();

    void mtKahyparTest();
#endif
    // analysis the different technology ratio
    double getAverageStdCellTechRatio(std::shared_ptr<PlacerBase> &pb_);

    double getAverageMacroTechRatio(std::shared_ptr<PlacerBase> &pb_);

    int getMacroStdcellAreaRatio(std::shared_ptr<PlacerBase> &pb_);

    double getDieUtilizeRatio(Die* die_);


private:
    // std::shared_ptr<Placer23b> placer23b_;
    // std::shared_ptr<Technology> technology_;
    // std::shared_ptr<PlacerBase> placerBase_;
    std::shared_ptr<Replace> replace_;
    // std::shared_ptr<> replace_;
  };
}

#endif