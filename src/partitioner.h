#ifndef __PARTITIONER__
#define __PARTITIONER__

#include "technology.h"
#include "placerBase.h"

namespace replace
{
  class Replace;

  class Partitioner
  {
  public:
    Partitioner(float targetDensity);
    ~Partitioner() = default;

    void partitioning(std::shared_ptr<PlacerBase> &pb_);
    void partitioning2(std::shared_ptr<PlacerBase> pb_);

    void hmetistest(std::shared_ptr<PlacerBase> pb_);

  private:
    // std::shared_ptr<Placer23b> placer23b_;
    // std::shared_ptr<Technology> technology_;
    // std::shared_ptr<PlacerBase> placerBase_;
    std::shared_ptr<Replace> replace_;
    // std::shared_ptr<> replace_;
  };
}

#endif
