#ifndef __REPLACE_LAYER_ASSIGNMENT__
#define __REPLACE_LAYER_ASSIGNMENT__

#include <memory>
#include <vector>
#include <cassert>

namespace replace{

    class PlacerBase;
    class NesterovBase;
    class GCell;

    class LayerAssignmenterVars
    {};

    class OverLap
    {
    public:
        OverLap(GCell* gcell1, GCell* gcell2);
        ~OverLap(){}

        void doOverLap();

        void init();
        void reset();

    private:
        std::vector<GCell*> gcells_;
        float overlapArea_;
    };

    class LayerAssignmenter
    {
    public:
        LayerAssignmenter();
        ~LayerAssignmenter();

        void setPlacerBase(const std::shared_ptr<PlacerBase>& pb) { pb_ = pb;}
        void setNesterovBase(const std::shared_ptr<NesterovBase>& nb) { nb_ = nb;}

        void doLayerAssignmenter();

        void init();
        void reset();

    private:
        std::shared_ptr<PlacerBase> pb_;
        std::shared_ptr<NesterovBase> nb_;

        std::vector<OverLap> overlapStor_;
        std::vector<GCell*> topdieGCells_;
        std::vector<GCell*> bottomdieGCells_;
    };

}

#endif