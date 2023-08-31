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
#if !defined(WIN32) && !defined(_WIN32)
    void mtPartitionInstance(std::shared_ptr<PlacerBase> &pb_);

    void partitionTest();

    void mtKahyparTest();
#endif
    // analysis the different technology
    double getAverageTechRatio(std::shared_ptr<PlacerBase> &pb_);

    int getMacroStdcellAreaRatio(std::shared_ptr<PlacerBase> &pb_);
  };
}

#endif
