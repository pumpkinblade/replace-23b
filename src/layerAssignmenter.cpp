#include "layerAssignmenter.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include "log.h"

namespace replace{
    OverLap::OverLap(GCell *gcell1, GCell *gcell2)
    {
        gcells_.push_back(gcell1);
        gcells_.push_back(gcell2);
    }
    LayerAssignmenter::LayerAssignmenter()
    {
    }
    LayerAssignmenter::~LayerAssignmenter()
    {
    }
    void LayerAssignmenter::doLayerAssignmenter()
    {
        // 初始化topdieGCells_和bottomdieGCells_
        topdieGCells_= std::vector<GCell*>();
        bottomdieGCells_ = std::vector<GCell*>();

        // 将同一个binGrid中的gcell进行随机划分到不同layer中
        for(BinGrid* bg : nb_->binGrids())
        {
            for (GCell *cell : bg->gCells())
            {
                // 随机划分到不同layer中
                int layer = rand() % 2;
                if (layer == 0){
                    topdieGCells_.emplace_back(cell);
                    LOG_INFO("topcell: ");
                }
                else{
                    bottomdieGCells_.emplace_back(cell);
                    LOG_INFO("bottomcell: ");
                }

            }
        }

        // 然后对每一个overlap进行layer assignment
    }
}