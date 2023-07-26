#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "placer23b.h"
#include "partitioner.h"
#include "replace.h"
#include <unordered_map>
#include <algorithm>

namespace replace
{
  Partitioner::Partitioner(float targetDensity)
  {
    std::unique_ptr<Replace> rp(new Replace(targetDensity));
    replace_ = std::move(rp);
  }

  void Partitioner::partitioning(std::shared_ptr<PlacerBase> &pb_)
  {
    // 先使用replace进行以此global placement
    replace_->setPlacerBase(pb_);
    replace_->doInitialPlace();
    replace_->doNesterovPlace("pregp");

    // TODO: move cell to bottom only where exists overlap
    LOG_TRACE("start partition");
    srand(10086);

    // create bottom nets
    std::vector<Net> bottomNets;
    const std::vector<Net *> &topNets = pb_->nets();
    bottomNets.reserve(pb_->nets().size());
    for (auto topnetptr : pb_->nets())
    {
      bottomNets.emplace_back();
    }
    std::unordered_map<Net *, Net *> topBotMap;
    for (int i = 0; i < bottomNets.size(); i++)
    {
      topBotMap.emplace(topNets[i], &bottomNets[i]);
    }

    // Short alias for top die and bottom die;
    Die &topdie = *pb_->die("top");
    Die &bottomdie = *pb_->die("bottom");
    std::vector<Instance*> macros;
    // 对replace布局后的结果进行layer assignment
    for (Instance *instance : pb_->insts())
    {
      // colloect all macros
      if(instance->isMacro()){
        macros.push_back(instance);
        continue;
      }
      float roll = (float)rand() / RAND_MAX;
      bool moveToBottom = roll >= 0.5 ? true : false;
      if (moveToBottom)
      {
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

        // TODO: Also, under different tech node, the pin offset from instance
        // is different
      }
    }

    std::sort(macros.begin(), macros.end(), [](const Instance* left,
    const Instance* right){ return left->size() > right->size(); });
    // A greedy macro partition method, which may be not optimal.
    int topMacroSize = 0;
    int bottomMacroSize = 0;
    for(Instance* macro : macros){
      if(topMacroSize > bottomMacroSize){
        topdie.removeInstance(macro);
        macro->setSize(*bottomdie.tech());
        bottomdie.addInstance(macro);

        for (Pin *pin : macro->pins())
        {
          auto topnet = pin->net();
          topnet->removePin(pin);
          topBotMap[topnet]->addPin(pin);
        }

        bottomMacroSize += macro->size();
      } else {
        topMacroSize += macro->size();
      }
    }
    LOG_INFO("partition macros, top: {} bottom: {}",
      topMacroSize, bottomMacroSize);


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
