#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "partitioner.h"
#include "replace.h"
#include <unordered_map>
#include <algorithm>
#include <unordered_set>


namespace replace
{
  Partitioner::Partitioner(float targetDensity)
  {
    std::unique_ptr<Replace> rp(new Replace(targetDensity));
    replace_ = std::move(rp);
  }

  // void Partitioner::partitioning(std::shared_ptr<PlacerBase> &pb_)
  // {
  //   // 先使用replace进行以此global placement
  //   replace_->setPlacerBase(pb_);
  //   replace_->doInitialPlace();
  //   replace_->doNesterovPlace("pregp");

  //   // TODO: move cell to bottom only where exists overlap
  //   LOG_TRACE("start partition");
  //   srand(10086);

  //   // create bottom nets
  //   std::vector<Net> bottomNets;
  //   const std::vector<Net *> &topNets = pb_->nets();
  //   bottomNets.reserve(pb_->nets().size());
  //   for (auto topnetptr : pb_->nets())
  //   {
  //     bottomNets.emplace_back();
  //   }
  //   std::unordered_map<Net *, Net *> topBotMap;
  //   for (int i = 0; i < bottomNets.size(); i++)
  //   {
  //     topBotMap.emplace(topNets[i], &bottomNets[i]);
  //   }

  //   // Short alias for top die and bottom die;
  //   Die &topdie = *pb_->die("top");
  //   Die &bottomdie = *pb_->die("bottom");
  //   std::vector<Instance*> macros;
  //   // 对replace布局后的结果进行layer assignment
  //   for (Instance *instance : pb_->insts())
  //   {
  //     // colloect all macros
  //     if(instance->isMacro()){
  //       macros.push_back(instance);
  //       continue;
  //     }
  //     float roll = (float)rand() / RAND_MAX;
  //     bool moveToBottom = roll >= 0.5 ? true : false;
  //     if (moveToBottom)
  //     {
  //       topdie.removeInstance(instance);

  //       // set instance size by bottom die technology
  //       instance->setSize(*bottomdie.tech());
  //       bottomdie.addInstance(instance);

  //       // When an instance is moved from top to bottom, nets connected
  //       // by the pins of this instance should be spilted to two
  //       // because current Net Class does not keep layer info of pin and
  //       // calculate HPWL assuming pins on the same layer.

  //       // So at the start, all nets are on the top layer, we name them
  //       // TOP NETS, and we create BOTTOM NETS, each bottom net
  //       // correponds to one top net. When we move instance to bottom,
  //       // we also remove pin in top net and add it to bottom net.
  //       for (Pin *pin : instance->pins())
  //       {
  //         auto topnet = pin->net();
  //         topnet->removePin(pin);
  //         topBotMap[topnet]->addPin(pin);
  //       }

  //       // TODO: Also, under different tech node, the pin offset from instance
  //       // is different
  //     }
  //   }

  //   std::sort(macros.begin(), macros.end(), [](const Instance* left,
  //   const Instance* right){ return left->size() > right->size(); });
  //   // A greedy macro partition method, which may be not optimal.
  //   int topMacroSize = 0;
  //   int bottomMacroSize = 0;
  //   for(Instance* macro : macros){
  //     if(topMacroSize > bottomMacroSize){
  //       topdie.removeInstance(macro);
  //       macro->setSize(*bottomdie.tech());
  //       bottomdie.addInstance(macro);

  //       for (Pin *pin : macro->pins())
  //       {
  //         auto topnet = pin->net();
  //         topnet->removePin(pin);
  //         topBotMap[topnet]->addPin(pin);
  //       }

  //       bottomMacroSize += macro->size();
  //     } else {
  //       topMacroSize += macro->size();
  //     }
  //   }
  //   LOG_DEBUG("partition macros, top: {} bottom: {}",
  //     topMacroSize, bottomMacroSize);


  //   // TODO: when this loop finish, we should clean empty top nets.
  //   //          and add terminal
  //   for (int i = 0; i < bottomNets.size(); i++)
  //   {
  //     auto topNet = topNets[i];
  //     Net &botNet = *topBotMap[topNet];
  //     // if bottom net is not empty
  //     if (botNet.pins().size() > 0)
  //     {
  //       // TODO: if top net is empty, remove this top net
  //       if (topNet->pins().size() == 0)
  //       {
  //         ;
  //       }
  //       else
  //       { // add terminal
  //         Instance &inst = pb_->emplaceInstance(false, false);
  //         inst.setSize(pb_->terminalSizeX(), pb_->terminalSizeY());
  //         inst.setLocation(2 * pb_->terminalSizeX(), 2 * pb_->terminalSizeY());
  //         inst.setFixed(false);
  //         // set inst's name using net name
  //         inst.setName(topNet->name());
  //         pb_->pinStor_.emplace_back();
  //         Pin &pinTop = pb_->pinStor_.emplace_back();
  //         Pin &pinBot = pb_->emplacePin();
  //         // add pin to net
  //         topNet->addPin(&pinTop);
  //         botNet.addPin(&pinBot);
  //         // add pin to instance
  //         inst.addPin(&pinTop);
  //         inst.addPin(&pinBot);
  //         // add instace to die
  //         pb_->die("terminal")->addInstance(&inst);
  //       }
  //       pb_->addNet(botNet);
  //     }
  //   }

  //   LOG_TRACE("finish partition");
  // }

  void moveDecide(int64_t topArea, int64_t botArea, int64_t topCap, int64_t botCap,
                  bool* moveToTop, bool* moveToBot)
  {
    bool topCanContain = topArea < topCap;
    bool botCanContain = botArea < botCap;

    *moveToTop = *moveToBot = false;
    if(topCanContain && botCanContain)
    {
      *moveToTop = ((float)rand() / RAND_MAX) < 0.5;
      *moveToBot = !(*moveToTop);
    }
    else if(topCanContain)
    {
      *moveToTop = true;
      *moveToBot = false;
    }
    else
    {
      *moveToBot = true;
      *moveToTop = false;
    }
  }

  void Partitioner::partitioning2(std::shared_ptr<PlacerBase> pb)
  {
    srand(114);

    Die* topdie = pb->die("top");
    Die* botdie = pb->die("bottom");

    // die capacity
    int64_t topCap = topdie->coreDx() * (long long)topdie->coreDy();
    topCap = static_cast<long long>(topCap * topdie->maxUtil());
    int64_t botCap = botdie->coreDx() * (long long)botdie->coreDy();
    botCap = static_cast<long long>(botCap * botdie->maxUtil());

    // cheat: using extId
    for(int i = 0; i < pb->insts().size(); i++)
    {
      pb->insts()[i]->setExtId(i);
    }
    std::vector<bool> hasAssigned(pb->insts().size(), false);
    std::vector<bool> isBot(pb->insts().size(), false);

    // using zjl's method to assign macros
    std::vector<Instance*> macros;
    for(Instance* inst : pb->insts())
    {
      if(inst->isMacro())
        macros.push_back(inst);
    }
    std::sort(macros.begin(), macros.end(), [](const Instance* left, const Instance* right)
              { return left->size() > right->size(); });
    // A greedy macro partition method, which may be not optimal.
    int topMacroSize = 0;
    int bottomMacroSize = 0;
    for(Instance* macro : macros)
    {
      if(topMacroSize > bottomMacroSize)
      {
        topdie->removeInstance(macro);
        macro->setSize(*botdie->tech());
        botdie->addInstance(macro);

        bottomMacroSize += macro->size();
        botCap -= macro->size();
        hasAssigned[macro->extId()] = true;
        isBot[macro->extId()] = true;
      }
      else
      {
        topMacroSize += macro->size();
        topCap -= macro->size();
        hasAssigned[macro->extId()] = true;
        isBot[macro->extId()] = false;
      }
    }
    LOG_DEBUG("partition macros, top: {} bottom: {}", topMacroSize, bottomMacroSize);

    // enumerating net
    for(Net* net : pb->nets())
    {
      std::unordered_set<Instance*> netInstSet;
      bool canPlaceTop = true;
      int64_t netAreaTop = 0;
      bool canPlaceBot = true;
      int64_t netAreaBot = 0;
      for(Pin* pin : net->pins())
      {
        Instance* inst = pin->instance();
        LibCell* libcell = botdie->tech()->libCell(inst->libCellId());
        canPlaceTop &= (!hasAssigned[inst->extId()]) || (!isBot[inst->extId()]);
        netAreaTop += hasAssigned[inst->extId()] ? 0 : (int64_t)inst->dx() * inst->dy();
        canPlaceBot &= (!hasAssigned[inst->extId()]) || (isBot[inst->extId()]);
        netAreaBot += hasAssigned[inst->extId()] ? 0 : (int64_t)libcell->sizeX() * libcell->sizeY();
      }

      bool toTop = false;
      bool toBot = false;
      if(canPlaceTop && netAreaTop < topCap && canPlaceBot && netAreaBot < botCap)
      {
        moveDecide(netAreaTop, netAreaBot, topCap, botCap, &toTop, &toBot);
      }
      else if(canPlaceTop && netAreaTop < topCap)
      {
        toTop = true;
      }
      else if(canPlaceBot && netAreaBot < botCap)
      {
        toBot = true;
      }

      if(toTop)
      {
        topCap -= netAreaTop;
        for(Pin* pin : net->pins())
        {
          Instance* inst = pin->instance();
          hasAssigned[inst->extId()] = true;
          isBot[inst->extId()] = false;
        }
      }
      else if(toBot)
      {
        botCap -= netAreaBot;
        for(Pin* pin : net->pins())
        {
          Instance* inst = pin->instance();
          if(hasAssigned[inst->extId()])
            continue;
          hasAssigned[inst->extId()] = true;
          isBot[inst->extId()] = true;
          inst->setSize(*botdie->tech());
          topdie->removeInstance(inst);
          botdie->addInstance(inst);
        }
      }
    }

    // for insts that hasn't been assigned
    for(Instance* inst : pb->insts())
    {
      if(hasAssigned[inst->extId()] == false)
      {
        auto instTopArea = inst->dx() * (long long)inst->dy();
        LibCell* botLibCell = botdie->tech()->libCell(inst->libCellId());
        auto instBotArea = botLibCell->sizeX() * (long long)botLibCell->sizeY();

        bool moveToTop, moveToBot;
        moveDecide(instTopArea, instBotArea, topCap, botCap, &moveToTop, &moveToBot);
        if(moveToBot)
        {
          botCap -= instBotArea;
          hasAssigned[inst->extId()] = true;
          isBot[inst->extId()] = true;

          inst->setSize(*botdie->tech());
          topdie->removeInstance(inst);
          botdie->addInstance(inst);
        }
        else
        {
          topCap -= instTopArea;
          isBot[inst->extId()] = false;
          hasAssigned[inst->extId()] = true;
        }
      }
    }

    // generate terminals
    int numTerms = 0;
    for(int i = 0, ie = pb->nets().size(); i < ie; i++)
    {
      Net* net = pb->nets()[i];
      bool hasTopPin = false;
      bool hasBotPin = false;
      for(Pin* pin : net->pins())
      {
        Instance* inst = pin->instance();
        hasTopPin |= !isBot[inst->extId()];
        hasBotPin |= isBot[inst->extId()];
      }
      numTerms += (hasTopPin && hasBotPin);
    }
    pb->extraInstStor_.reserve(numTerms);
    pb->extraNetStor_.reserve(numTerms);
    pb->extraPinStor_.reserve(2*numTerms);
    pb->insts_.reserve(pb->insts_.size() + numTerms);
    pb->nets_.reserve(pb->nets_.size() + numTerms);
    pb->pins_.reserve(pb->pins_.size() + 2 * numTerms);

    for(int i = 0, ie = pb->nets().size(); i < ie; i++)
    {
      Net* net = pb->nets()[i];
      bool hasTopPin = false;
      bool hasBotPin = false;
      for(Pin* pin : net->pins())
      {
        Instance* inst = pin->instance();
        hasTopPin |= !isBot[inst->extId()];
        hasBotPin |= isBot[inst->extId()];
      }
      // need to split net and generate net
      if(hasTopPin && hasBotPin)
      {
        // generate terminal
        pb->extraInstStor_.emplace_back();
        Instance* term = &pb->extraInstStor_.back();
        term->setFixed(false);
        term->setMacro(false);
        term->setSize(pb->terminalSizeX(), pb->terminalSizeY());
        term->setLocation(2 * pb->terminalSizeX(), 2 * pb->terminalSizeY());
        term->setName(pb->netNameStor_[i]);
        pb->die("terminal")->addInstance(term);
        pb->insts_.push_back(term);

        // generate pin for term
        pb->extraPinStor_.emplace_back();
        Pin* pinTop = &pb->extraPinStor_.back();
        pb->extraPinStor_.emplace_back();
        Pin* pinBot = &pb->extraPinStor_.back();
        term->addPin(pinTop);
        term->addPin(pinBot);
        pb->pins_.push_back(pinTop);
        pb->pins_.push_back(pinBot);
        // add instance to die
        
        // generate another net
        pb->extraNetStor_.emplace_back();
        Net* net2 = &pb->extraNetStor_.back();
        pb->nets_.push_back(net2);
        // add pin to net
        for(Pin* pin : net->pins())
        {
          int id = pin->instance()->extId();
          if(isBot[id])
          {
            net->removePin(pin);
            net2->addPin(pin);
          }
        }
        net->addPin(pinTop);
        net2->addPin(pinBot);
      }
    }

    // clear net name
    pb->netNameStor_ = std::vector<string>();
  }
}
