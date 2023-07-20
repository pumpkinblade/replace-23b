#ifndef __PARTITIONER__
#define __PARTITIONER__

#include "technology.h"
#include "placerBase.h"
#include "placer23b.h"

namespace replace{
    class Replace;

    class P23bToBaseConverter{
        public:
            P23bToBaseConverter() = default;
            ~P23bToBaseConverter() = default;

            std::shared_ptr<PlacerBase> placer23bToPlaceBase(std::shared_ptr<Placer23b> placer23b_);

    };



    class Partitioner{
        public:
            Partitioner(float targetDensity);
            void partitioning(std::shared_ptr<PlacerBase>& pb_);
            ~Partitioner() = default;

            // 将instance23b转换为libcell
            std::vector<LibCell*> InstanceToLibCells(std::vector<Instance23b*> instances, Technology* technology);
            // 将Placer23b根据partitioning的要求转换为PlacerBase
            PlacerBase* partitionToPlacerBase(Placer23b* placer23b);

        private:
            // std::shared_ptr<Placer23b> placer23b_;
            // std::shared_ptr<Technology> technology_;
            // std::shared_ptr<PlacerBase> placerBase_;
            std::shared_ptr<Replace> replace_;
            // std::shared_ptr<> replace_;


    };
}

#endif
