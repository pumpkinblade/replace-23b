#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "placer23b.h"
#include "partitioner.h"
#include "replace.h"

namespace replace{

    std::shared_ptr<PlacerBase> P23bToBaseConverter::placer23bToPlaceBase(std::shared_ptr<Placer23b> placer23b_){
        // auto pb = std::make_shared<PlacerBase>();
        // Technology* topDieTechnology = placer23b_->topDieTechnology();
        // // Process die and rows
        // pb->dieStor_.emplace_back(placer23b_->topDie());
        // pb->dies_.push_back(&pb->dieStor_.back());
        
        // // 给placebase预留空间
        // LOG_TRACE("Process def component");
        // std::unordered_map<std::string, int> instExtIds;
        // pb->instStor_.reserve(placer23b_->insts().size());
        // for (auto &inst : placer23b_->insts())
        // {
        //     // 在Technology中找到对应的LibCell
        //     auto libCell = topDieTechnology->cell(inst->libCellName());
        //     // 将LibCell转换为Instance
        //     Instance instance;
        //     instance.setSize(libCell->sizeX(), libCell->sizeY());
        //     bool isMacro = libCell->isMacro();
        //     instance.setMacro(isMacro);
        //     instance.setFixed(0);
        //     // 将Instance加入到PlacerBase中
        //     pb->instStor_.push_back(instance);
            
        // }
        // LOG_TRACE("Process def net");
        // // Process def net
        // pb->netStor_.reserve(placer23b_->nets().size());
        // for (auto &net : placer23b_->nets()){
        //     pb->netStor_.emplace_back();
        //     for(auto &pin23b : net->pins()){
        //         Pin pin;
        //         // pin.setInstance(pin23b->instName());
        //     }
        // }

    }

    Partitioner::Partitioner(float targetDensity){
        std::unique_ptr<Replace> rp(new Replace(targetDensity));
        replace_ = std::move(rp);
    }


    void Partitioner::partitioning(std::shared_ptr<PlacerBase>& pb_){
        // 先使用replace进行以此global placement
        replace_->setPlacerBase(pb_);
        replace_->doInitialPlace();
        replace_->doNesterovPlace();

        // TODO: move cell to bottom only where exists overlap
        LOG_TRACE("start partition");
        srand(10086);
        // 对replace布局后的结果进行layer assignment
        for (auto instance : pb_->insts()){
            float roll = (float) rand()/RAND_MAX;
            bool moveToBottom = roll >= 0.5 ? true : false;
            if(moveToBottom){
                pb_->die("top")->removeInstance(instance);
                pb_->die("bottom")->addInstance(instance);
            }
        }

        // TODO: add terminal instance to terminal die


        LOG_TRACE("finish partition");
    }
}
