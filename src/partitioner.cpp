#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "placer23b.h"
#include "partitioner.h"
#include "replace.h"
#include <unordered_map>


namespace replace
{

  std::shared_ptr<PlacerBase> P23bToBaseConverter::placer23bToPlaceBase(std::shared_ptr<Placer23b> placer23b_)
  {
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
    return nullptr;
  }

  Partitioner::Partitioner(float targetDensity)
  {
    std::unique_ptr<Replace> rp(new Replace(targetDensity));
    replace_ = std::move(rp);
  }

  void Partitioner::partitioning(std::shared_ptr<PlacerBase> &pb_)
  {
    assert(&pb_->instStor_[0] == pb_->insts_[0]);
    // 先使用replace进行以此global placement
    replace_->setPlacerBase(pb_);
    replace_->doInitialPlace();
    replace_->doNesterovPlace("pregp");
    assert(&pb_->instStor_[0] == pb_->insts_[0]);

    // TODO: move cell to bottom only where exists overlap
    LOG_TRACE("start partition");
    srand(10086);

    // create bottom nets
    std::vector<Net> bottomNets;
    const std::vector<Net *> &topNets = pb_->nets();
    bottomNets.reserve(pb_->nets().size());
    assert(&pb_->instStor_[0] == pb_->insts_[0]);
    for (auto topnetptr : pb_->nets())
    {
      bottomNets.emplace_back();
    }
    std::unordered_map<Net *, Net *> topBotMap;
    for (int i = 0; i < bottomNets.size(); i++)
    {
      topBotMap.emplace(topNets[i], &bottomNets[i]);
    }

    // 对replace布局后的结果进行layer assignment
    for (Instance *instance : pb_->insts())
    {
      assert(&pb_->instStor_[0] == pb_->insts_[0]);
      float roll = (float)rand() / RAND_MAX;
      bool moveToBottom = roll >= 0.5 ? true : false;
      if (moveToBottom)
      {
        Die &topdie = *pb_->die("top");
        Die &bottomdie = *pb_->die("bottom");

        topdie.removeInstance(instance);

        // set instance size by bottom die technology
        instance->setSize(*bottomdie.tech());
        bottomdie.addInstance(instance);

        // When an instance is moved from top to bottom, nets connected
        // by the pins of this instance should be spilted to two
        // because current Net Class does not keep layer info of pin and
        // calculate HPWL assuming pins on the same layer.

        // So at the start, all nets are on the top layer, we name them
        // TOP NETS, and we create BOTTOM NETS, each bottom net
        // correponds to one top net. When we move instance to bottom,
        // we also remove pin in top net and add it to bottom net.
        for (Pin *pin : instance->pins())
        {
          auto topnet = pin->net();
          topnet->removePin(pin);
          topBotMap[topnet]->addPin(pin);
        }

        // Also, under different tech node, the pin offset from instance
        // is different
      }
    }

    // TODO: when this loop finish, we should clean empty top nets.
    //          and add terminal
    for (int i = 0; i < bottomNets.size(); i++)
    {
      auto topNet = topNets[i];
      Net &botNet = *topBotMap[topNet];
      // if bottom net is not empty
      if (botNet.pins().size() > 0)
      {
        // TODO: if top net is empty, remove this top net
        if (topNet->pins().size() == 0)
        {
          ;
        }
        else
        { // add terminal
          Instance &inst = pb_->emplaceInstance(false, false);
          inst.setSize(pb_->terminalSizeX(), pb_->terminalSizeY());
          inst.setLocation(2 * pb_->terminalSizeX(), 2 * pb_->terminalSizeY());
          inst.setFixed(false);
          // set inst's name using net name
          inst.setName(topNet->name());
          Pin &pinTop = pb_->emplacePin();
          Pin &pinBot = pb_->emplacePin();
          // add pin to net
          topNet->addPin(&pinTop);
          botNet.addPin(&pinBot);
          // add pin to instance
          inst.addPin(&pinTop);
          inst.addPin(&pinBot);
          // add instace to die
          pb_->die("terminal")->addInstance(&inst);
        }
        pb_->addNet(botNet);
      }
    }

    LOG_TRACE("finish partition");
  }
}
