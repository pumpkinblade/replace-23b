#ifndef __RANDOM_PARTI__
#define __RANDOM_PARTI__

#include "placer23b.h"
#include "placerBase.h"

namespace replace{
class RandomPartitioner{
public:
    RandomPartitioner(Placer23b& p23): p23_(p23){};
    PlacerBase doPartition();
private:
    Placer23b& p23_;
};
}

#endif
