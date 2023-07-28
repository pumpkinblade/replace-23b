#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "placer23b.h"
#include "partitioner.h"
#include "replace.h"
#include <unordered_map>
#include <unordered_set>


namespace replace
{
  Partitioner::Partitioner(float targetDensity)
  {
    std::unique_ptr<Replace> rp(new Replace(targetDensity));
    replace_ = std::move(rp);
  }

  void Partitioner::partitioning(std::shared_ptr<PlacerBase> &pb_)
  {
    // assert(&pb_->instStor_[0] == pb_->insts_[0]);
    // // 先使用replace进行以此global placement
    // replace_->setPlacerBase(pb_);
    // replace_->doInitialPlace();
    // replace_->doNesterovPlace("pregp");
    // assert(&pb_->instStor_[0] == pb_->insts_[0]);

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
    srand(114514);

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

    // for insts that hasn't been assigned
    for(Instance* inst : pb->insts())
    {
      if(hasAssigned[inst->extId()] == false)
      {
        auto instTopArea = inst->dx() * (long long)inst->dy();
        LibCell* botLibCell = botdie->tech()->libCell(inst->libCellName());
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
    for(Net* net : pb->nets())
    {
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
        // add terminal
        Instance &term = pb->emplaceInstance(false, false);
        term.setSize(pb->terminalSizeX(), pb->terminalSizeY());
        term.setLocation(2 * pb->terminalSizeX(), 2 * pb->terminalSizeY());
        term.setFixed(false);
        // set term's name using net name
        term.setName(net->name());
        // generate pin for term
        Pin &pinTop = pb->emplacePin();
        Pin &pinBot = pb->emplacePin();
        term.addPin(&pinTop);
        term.addPin(&pinBot);
        // add instance to die
        pb->instNameMap_.emplace(term.name(), &term);
        pb->die("terminal")->addInstance(&term);
        
        // generate another net
        pb->netStor_.emplace_back();
        Net* net2 = &pb->netStor_.back();
        pb->netNameMap_.emplace(net->name() + "2", net2);
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
        net->addPin(&pinTop);
        net2->addPin(&pinBot);
      }
    }
  }
}
