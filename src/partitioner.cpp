#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "partitioner.h"
#include "replace.h"
#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include <libkahypar.h>
#include <iostream>
#include <vector>
#include <libmtkahypar.h>
#include <thread>
#include <memory>
#include <sys/time.h>
#include <iomanip>

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
      float thres = static_cast<float>(botArea) / static_cast<float>(botArea + topArea);
      *moveToTop = ((float)rand() / RAND_MAX) < thres;
      *moveToBot = !(*moveToTop);
    }
    else if(topCanContain)
    {
      *moveToBot = false;
      *moveToTop = true;
    }
    else
    {
      *moveToBot = true;
      *moveToTop = false;
    }
  }

  inline int64_t instCap(int dx, int dy, bool isMacro, Die* die)
  {
    if(isMacro)
      return dx * (dy + die->rowHeight() - 1) / die->rowHeight() * die->rowHeight();
    else
      return dx * dy;
  }

  void Partitioner::partitioning2(std::shared_ptr<PlacerBase> pb)
  {
    srand(514);

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
        canPlaceBot &= (!hasAssigned[inst->extId()]) || (isBot[inst->extId()]);
        netAreaTop += hasAssigned[inst->extId()] ? 0 : instCap(inst->dx(), inst->dy(), inst->isMacro(), topdie);
        netAreaBot += hasAssigned[inst->extId()] ? 0 : instCap(libcell->sizeX(), libcell->sizeY(), inst->isMacro(), botdie);
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
    LOG_INFO("TopDie remaining capactiy: {}", topCap);
    LOG_INFO("BotDie remaining capactiy: {}", botCap);
    LOG_INFO("numTerms : {}", numTerms);
    if (topCap < 0 || botCap < 0)
    {
      LOG_CRITICAL("Violated max utilization!!");
      exit(0);
    }
    pb->extraInstStor_.reserve(numTerms);
    pb->extraNetStor_.reserve(numTerms);
    pb->extraPinStor_.reserve(2 * numTerms);
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


  void Partitioner::partitionInstance(std::shared_ptr<PlacerBase> &pb_){
    // convert the instance and net to hypergraph
    LOG_INFO("start partition");
    // cheat: using extId
    for(int i = 0; i < pb_->insts().size(); i++)
    {
      pb_->insts()[i]->setExtId(i);
    }
    int instanceNum=pb_->insts_.size();
    int netNum=pb_->nets().size();
    const kahypar_hypernode_id_t num_vertices = instanceNum;
    const kahypar_hyperedge_id_t num_hyperedges = netNum;

    std::unique_ptr<kahypar_hyperedge_weight_t[]> hyperedge_weights = std::make_unique<kahypar_hyperedge_weight_t[]>(netNum);

    for(int i=0;i<netNum;i++){
      hyperedge_weights[i] = 1;
    }

    // calculate the length of hyperedges
    int hyperedgesLength = 0;
    for(Net* curNet : pb_->nets()){
      int netLength = curNet->pins().size();
      hyperedgesLength += netLength;
    }

    // initial hyperedges
    std::unique_ptr<kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<kahypar_hyperedge_id_t[]>(hyperedgesLength);
    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(netNum+1);
    int indicesIdx = 0;
    for(int i=0;i<netNum;i++){
      hyperedge_indices[i]=indicesIdx;
      Net* curNet = pb_->nets()[i];
      for(int j=0 ; j<curNet->pins().size(); j++){
        Pin* curPin = curNet->pins()[j];
        int nodeIndex=curPin->instance()->extId();
        // std::cout<<"extid: "<<nodeIndex<<std::endl;
        hyperedges[indicesIdx]=nodeIndex;
        indicesIdx++;
      }
    }
    hyperedge_indices[netNum]=indicesIdx-1;
    const double imbalance = 0.1;
    const kahypar_partition_id_t k = 2;

    kahypar_hyperedge_weight_t objective = 0;

    std::vector<kahypar_partition_id_t> partition(num_vertices, -1);

    kahypar_context_t* context = kahypar_context_new();
    kahypar_configure_context_from_file(context, "/home/workspace/replace-parker/replace-23b/test/cut_kKaHyPar_sea20.ini");
    
    kahypar_set_seed(context, 42);

    kahypar_partition(num_vertices, num_hyperedges,
                      imbalance, k,
                      /*vertex_weights */ nullptr, hyperedge_weights.get(),
                      hyperedge_indices.get(), hyperedges.get(),
                      &objective, context, partition.data());

    for(int i = 0; i != num_vertices; ++i) {
      std::cout << i << ":" << partition[i] << std::endl;
    }

    // move the instance to the bottom die
    // create bottom nets
    std::cout<<"---------------------1"<<std::endl;

    // Short alias for top die and bottom die;
    Die* topdie = pb_->die("top");
    Die* botdie = pb_->die("bottom");
    std::vector<bool> isBot(pb_->insts().size(), false);
    std::vector<Instance*> macros;
    for(Instance* inst : pb_->insts())
    {
      // collect all macro
      if(inst->isMacro()){
        macros.push_back(inst);
        continue;
      }
      auto instTopArea = inst->dx() * (long long)inst->dy();
      LibCell* botLibCell = botdie->tech()->libCell(inst->libCellId());
      auto instBotArea = botLibCell->sizeX() * (long long)botLibCell->sizeY();

      bool moveToTop, moveToBot;
      auto extid=inst->extId();
      if (partition[extid] == 0){
        moveToTop =  true;
      }else{
        moveToTop = false;
        moveToBot = true;
      }
      if(moveToBot)
      {
        isBot[inst->extId()] = true;
        inst->setSize(*botdie->tech());
        topdie->removeInstance(inst);
        botdie->addInstance(inst);
      }
      else
      {
        isBot[inst->extId()] = false;
      }
    }
    kahypar_context_free(context);
    std::cout<<"---------------------2"<<std::endl;
    std::sort(macros.begin(), macros.end(), [](const Instance* left,
    const Instance* right){ return left->size() > right->size(); });
    // A greedy macro partition method, which may be not optimal.
    int topMacroSize = 0;
    int bottomMacroSize = 0;
    for(Instance* macro : macros){
      if(topMacroSize > bottomMacroSize){
        isBot[macro->extId()] = true;
        macro->setSize(*botdie->tech());
        topdie->removeInstance(macro);
        botdie->addInstance(macro);
        bottomMacroSize += macro->size();
      }else {
        topMacroSize += macro->size();
      }
    }
    std::cout<<"---------------------3"<<std::endl;
    // generate terminals
    int numTerms = 0;
    for(int i = 0, ie = pb_->nets().size(); i < ie; i++)
    {
      Net* net = pb_->nets()[i];
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
    pb_->extraInstStor_.reserve(numTerms);
    pb_->extraNetStor_.reserve(numTerms);
    pb_->extraPinStor_.reserve(2 * numTerms);
    pb_->insts_.reserve(pb_->insts_.size() + numTerms);
    pb_->nets_.reserve(pb_->nets_.size() + numTerms);
    pb_->pins_.reserve(pb_->pins_.size() + 2 * numTerms);

    for(int i = 0, ie = pb_->nets().size(); i < ie; i++)
    {
      Net* net = pb_->nets()[i];
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
        pb_->extraInstStor_.emplace_back();
        Instance* term = &pb_->extraInstStor_.back();
        term->setFixed(false);
        term->setMacro(false);
        term->setSize(pb_->terminalSizeX(), pb_->terminalSizeY());
        term->setLocation(2 * pb_->terminalSizeX(), 2 * pb_->terminalSizeY());
        term->setName(pb_->netNameStor_[i]);
        pb_->die("terminal")->addInstance(term);
        pb_->insts_.push_back(term);

        // generate pin for term
        pb_->extraPinStor_.emplace_back();
        Pin* pinTop = &pb_->extraPinStor_.back();
        pb_->extraPinStor_.emplace_back();
        Pin* pinBot = &pb_->extraPinStor_.back();
        term->addPin(pinTop);
        term->addPin(pinBot);
        pb_->pins_.push_back(pinTop);
        pb_->pins_.push_back(pinBot);
        // add instance to die
        
        // generate another net
        pb_->extraNetStor_.emplace_back();
        Net* net2 = &pb_->extraNetStor_.back();
        pb_->nets_.push_back(net2);
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
    pb_->netNameStor_ = std::vector<string>();
    std::cout<<"---------------------4"<<std::endl;
    
    LOG_INFO("finish partition");
  }

  void Partitioner::mtPartitionInstance(std::shared_ptr<PlacerBase> &pb_){
    // convert the instance and net to hypergraph
    LOG_INFO("start partition");
    // cheat: using extId
    for(int i = 0; i < pb_->insts().size(); i++)
    {
      pb_->insts()[i]->setExtId(i);
    }
    int instanceNum=pb_->insts_.size();
    int netNum=pb_->nets().size();
    const mt_kahypar_hypernode_id_t num_vertices = instanceNum;
    const mt_kahypar_hyperedge_id_t num_hyperedges = netNum;

    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> hyperedge_weights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(netNum);

    for(int i=0;i<netNum;i++){
      hyperedge_weights[i] = 1;
    }
    // get macro and stdcell ratio
    // int macroStdcellAreaRatio=this->getMacroStdcellAreaRatio(pb_);
    long vertex_weights_sum=0;
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> vertex_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(instanceNum);
    for(int i=0;i<instanceNum;i++){
      // Instance* inst_=pb_->insts()[i];
      // if(inst_->isMacro()){
      //   vertex_weights[i] = macroStdcellAreaRatio;
      //   vertex_weights_sum+=macroStdcellAreaRatio;
      // }
      // else{
      //   vertex_weights[i] = 1;
      //   vertex_weights_sum+=1;
      // }
      vertex_weights[i] = 1;
      vertex_weights_sum+=1;
    }

    // calculate the length of hyperedges
    int hyperedgesLength = 0;
    for(Net* curNet : pb_->nets()){
      int netLength = curNet->pins().size();
      hyperedgesLength += netLength;
    }

    // initial hyperedges
    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(hyperedgesLength);
    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(netNum+1);
    int indicesIdx = 0;
    for(int i=0;i<netNum;i++){
      hyperedge_indices[i]=indicesIdx;
      Net* curNet = pb_->nets()[i];
      for(int j=0 ; j<curNet->pins().size(); j++){
        Pin* curPin = curNet->pins()[j];
        int nodeIndex=curPin->instance()->extId();
        // std::cout<<"extid: "<<nodeIndex<<std::endl;
        hyperedges[indicesIdx]=nodeIndex;
        indicesIdx++;
      }
    }
    hyperedge_indices[netNum]=indicesIdx-1;
    const mt_kahypar_partition_id_t k = 2;

    mt_kahypar_hyperedge_weight_t objective = 0;

    mt_kahypar_context_t* context = mt_kahypar_context_new();
    mt_kahypar_load_preset(context, DEFAULT);
    
    const mt_kahypar_partition_id_t num_blocks=2;
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> init_block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
    
    // get average technology ratio
    double averageRatio=this->getAverageTechRatio(pb_);
    if(averageRatio==1){
      double ratioA=0.5;
      std::cout<<"averageRatio: "<<averageRatio<<std::endl;
      init_block_weights[0]=int(ratioA*vertex_weights_sum*1.1);
      init_block_weights[1]=int(ratioA*vertex_weights_sum*1.1);
    }
    else{
      double ratioA=averageRatio/(averageRatio+1);
      double ratioB=1/(averageRatio+1);
      string topTech=pb_->die("top")->tech()->name();
      string botTech=pb_->die("bottom")->tech()->name();
      string Tech1Name=pb_->techs()[0]->name();
      string Tech2Name=pb_->techs()[1]->name();
      if(topTech==Tech1Name){
        init_block_weights[1]=int(ratioA*vertex_weights_sum*1.1);
        init_block_weights[0]=int(ratioB*vertex_weights_sum*1.1);
      }else if(topTech==Tech2Name){
        init_block_weights[0]=int(ratioA*vertex_weights_sum*1.1);
        init_block_weights[1]=int(ratioB*vertex_weights_sum*1.1);
      }
    }

    std::cout<<"init_block_weights[0]: "<<init_block_weights[0]<<std::endl;
    std::cout<<"init_block_weights[1]: "<<init_block_weights[1]<<std::endl;

    mt_kahypar_set_individual_target_block_weights(context, num_blocks, init_block_weights.get());
    mt_kahypar_set_partitioning_parameters(context,
      num_blocks /* number of blocks */, 0.1 /* imbalance parameter */,
      CUT /* objective function */);
    mt_kahypar_set_seed(42 /* seed */);

    mt_kahypar_set_context_parameter(context, VERBOSE, "1");

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(DEFAULT, num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), hyperedge_weights.get(), vertex_weights.get());
    // mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(DEFAULT, num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), hyperedge_weights.get(), nullptr);
    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);

    // Partition Hypergraph
    mt_kahypar_partitioned_hypergraph_t partitioned_hg =
      mt_kahypar_partition(hypergraph, context);
    
    // Extract Partition
    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
    mt_kahypar_get_partition(partitioned_hg, partition.get());

    // Extract Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
    mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

    // Compute Metrics
    const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
    const double km1 = mt_kahypar_km1(partitioned_hg);

    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds*1e-6;

    // move the instance to the bottom die
    // create bottom nets
    LOG_INFO("---------------------");
    std::cout<<"---------------------1"<<std::endl;

    // Output Results
    std::cout << "Partitioning Results:" << std::endl;
    std::cout << "Imbalance         = " << imbalance << std::endl;
    std::cout << "Km1               = " << km1 << std::endl;
    std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
    std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;
    std::cout << "Run time = " << std::setprecision(3)<< elapsed <<"s"<< std::endl;

    // Short alias for top die and bottom die;
    Die* topdie = pb_->die("top");
    Die* botdie = pb_->die("bottom");
    std::vector<bool> isBot(pb_->insts().size(), false);
    std::vector<Instance*> macros;
    for(Instance* inst : pb_->insts())
    {
      // collect all macro
      if(inst->isMacro()){
        macros.push_back(inst);
        continue;
      }
      // auto instTopArea = inst->dx() * (long long)inst->dy();
      // LibCell* botLibCell = botdie->tech()->libCell(inst->libCellId());
      // auto instBotArea = botLibCell->sizeX() * (long long)botLibCell->sizeY();

      auto extid=inst->extId();
      if (partition[extid] == 1){
        isBot[extid] = true;
        inst->setSize(*botdie->tech());
        topdie->removeInstance(inst);
        botdie->addInstance(inst);
      }else{
        isBot[extid] = false;
      }

    }
    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
    std::cout<<"---------------------2"<<std::endl;
    std::sort(macros.begin(), macros.end(), [](const Instance* left,
    const Instance* right){ return left->size() > right->size(); });
    // A greedy macro partition method, which may be not optimal.
    int topMacroSize = 0;
    int bottomMacroSize = 0;
    for(Instance* macro : macros){
      if(topMacroSize > bottomMacroSize){
        isBot[macro->extId()] = true;
        macro->setSize(*botdie->tech());
        topdie->removeInstance(macro);
        botdie->addInstance(macro);
        bottomMacroSize += macro->size();
      }else {
        topMacroSize += macro->size();
      }
    }
    std::cout<<"---------------------3"<<std::endl;
    // generate terminals
    int numTerms = 0;
    for(int i = 0, ie = pb_->nets().size(); i < ie; i++)
    {
      Net* net = pb_->nets()[i];
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
    pb_->extraInstStor_.reserve(numTerms);
    pb_->extraNetStor_.reserve(numTerms);
    pb_->extraPinStor_.reserve(2 * numTerms);
    pb_->insts_.reserve(pb_->insts_.size() + numTerms);
    pb_->nets_.reserve(pb_->nets_.size() + numTerms);
    pb_->pins_.reserve(pb_->pins_.size() + 2 * numTerms);

    for(int i = 0, ie = pb_->nets().size(); i < ie; i++)
    {
      Net* net = pb_->nets()[i];
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
        pb_->extraInstStor_.emplace_back();
        Instance* term = &pb_->extraInstStor_.back();
        term->setFixed(false);
        term->setMacro(false);
        term->setSize(pb_->terminalSizeX(), pb_->terminalSizeY());
        term->setLocation(2 * pb_->terminalSizeX(), 2 * pb_->terminalSizeY());
        term->setName(pb_->netNameStor_[i]);
        pb_->die("terminal")->addInstance(term);
        pb_->insts_.push_back(term);

        // generate pin for term
        pb_->extraPinStor_.emplace_back();
        Pin* pinTop = &pb_->extraPinStor_.back();
        pb_->extraPinStor_.emplace_back();
        Pin* pinBot = &pb_->extraPinStor_.back();
        term->addPin(pinTop);
        term->addPin(pinBot);
        pb_->pins_.push_back(pinTop);
        pb_->pins_.push_back(pinBot);
        // add instance to die
        
        // generate another net
        pb_->extraNetStor_.emplace_back();
        Net* net2 = &pb_->extraNetStor_.back();
        pb_->nets_.push_back(net2);
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
    pb_->netNameStor_ = std::vector<string>();
    std::cout<<"---------------------4"<<std::endl;
    
    LOG_INFO("finish partition");
} 
  
  void Partitioner::mtKahyparTest(){
    // Initialize thread pool
    // mt_kahypar_initialize_thread_pool(
    //   std::thread::hardware_concurrency() /* use all available cores */,
    //   true /* activate interleaved NUMA allocation policy */ );

    mt_kahypar_initialize_thread_pool(
      8 /* use all available cores */,
      true /* activate interleaved NUMA allocation policy */ );

    // Setup partitioning context
    mt_kahypar_context_t* context = mt_kahypar_context_new();
    mt_kahypar_load_preset(context, DEFAULT /* corresponds to MT-KaHyPar-D */);
    // In the following, we partition a hypergraph into two blocks
    // with an allowed imbalance of 3% and optimize the connective metric (KM1)
    mt_kahypar_set_partitioning_parameters(context,
      2 /* number of blocks */, 0.03 /* imbalance parameter */,
      KM1 /* objective function */);
    mt_kahypar_set_seed(42 /* seed */);
    // Enable logging
    mt_kahypar_set_context_parameter(context, VERBOSE, "1");

    // Load Hypergraph for DEFAULT preset
    // mt_kahypar_hypergraph_t hypergraph =
    //   mt_kahypar_read_hypergraph_from_file(
    //     "path/to/hypergraph/file", DEFAULT, HMETIS /* file format */);

    // Create Hypergraph
    const mt_kahypar_hypernode_id_t num_vertices = 5000;
    const mt_kahypar_hyperedge_id_t num_hyperedges = 1000;
    const int maxNodesPerHyperedge = 400;
    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> hyperedge_weights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(num_hyperedges);

    for(int i=0; i<num_hyperedges; i++){
      hyperedge_weights[i] = 1;
    }
    std::vector<int> hyperedge_weights_vec;
    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(num_hyperedges+1);
    srand(time(0));
    int indicesIdx = 0;
    for(int i=0; i<num_hyperedges; i++){
      int numNodesInHyperedge = rand() % maxNodesPerHyperedge + 1;
      hyperedge_indices[i]=indicesIdx;
      std::vector<int> hyperedge_edge_tmp;
      for(int j=0; j<numNodesInHyperedge; j++){
        int nodesIdxInHyperedge = rand() % num_vertices;
        hyperedge_edge_tmp.push_back(nodesIdxInHyperedge);
        indicesIdx++;
      }
      sort(hyperedge_edge_tmp.begin(),hyperedge_edge_tmp.end());
      hyperedge_weights_vec.insert(hyperedge_weights_vec.end(),hyperedge_edge_tmp.begin(),hyperedge_edge_tmp.end());
    }
    hyperedge_indices[num_hyperedges]=indicesIdx-1;
    int hyperedge_len=hyperedge_weights_vec.size();

    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<mt_kahypar_hyperedge_id_t[]>(hyperedge_len);

    for(int i=0;i<hyperedge_len;i++){
      hyperedges[i]=hyperedge_weights_vec[i];
    }

    mt_kahypar_hypergraph_t hypergraph = mt_kahypar_create_hypergraph(DEFAULT, num_vertices, num_hyperedges, hyperedge_indices.get(), hyperedges.get(), hyperedge_weights.get(), nullptr);

    // Start measuring time
    struct timeval begin, end;
    gettimeofday(&begin, 0);

    // Partition Hypergraph
    mt_kahypar_partitioned_hypergraph_t partitioned_hg =
      mt_kahypar_partition(hypergraph, context);

    // Extract Partition
    std::unique_ptr<mt_kahypar_partition_id_t[]> partition =
      std::make_unique<mt_kahypar_partition_id_t[]>(mt_kahypar_num_hypernodes(hypergraph));
    mt_kahypar_get_partition(partitioned_hg, partition.get());

    // Extract Block Weights
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> block_weights =
      std::make_unique<mt_kahypar_hypernode_weight_t[]>(2);
    mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

    // Compute Metrics
    const double imbalance = mt_kahypar_imbalance(partitioned_hg, context);
    const double km1 = mt_kahypar_km1(partitioned_hg);

    // Stop measuring time and calculate the elapsed time
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds*1e-6;

    // Output Results
    std::cout << "Partitioning Results:" << std::endl;
    std::cout << "Imbalance         = " << imbalance << std::endl;
    std::cout << "Km1               = " << km1 << std::endl;
    std::cout << "Weight of Block 0 = " << block_weights[0] << std::endl;
    std::cout << "Weight of Block 1 = " << block_weights[1] << std::endl;
    std::cout << "Run time = " << std::setprecision(3)<< elapsed <<"s"<< std::endl;

    mt_kahypar_free_context(context);
    mt_kahypar_free_hypergraph(hypergraph);
    mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  }

  double Partitioner::getAverageTechRatio(std::shared_ptr<PlacerBase> &pb_)
  {
    // get top die tech
    // Technology* ta=pb_->dies()[0]->tech();
    int technum = pb_->techs().size();
    if(technum==1){
      return double(1);
    }
    Technology* ta = pb_->techs()[0];
    // get bottom die tech
    // Technology* tb=pb_->dies()[2]->tech();
    Technology* tb = pb_->techs()[1];
    // traverse libcell in each tech to get the ratio of them
    int sizea=ta->libCells().size();
    int sizeb=tb->libCells().size();
    std::vector<double> ratiolist(sizea);
    for(int i=0;i<sizea;i++){
      LibCell* libcella=ta->libCells()[i];
      int ida=libcella->id();
      LibCell* libcellb=tb->libCell(ida);
      int areaA=libcella->sizeX()*libcella->sizeY();
      int areaB=libcellb->sizeX()*libcellb->sizeY();
      double tmpRatio=(double)areaA/(double)areaB;
      ratiolist[i]=tmpRatio;
    }
    // calculate average ratio
    double averageRatio=ratiolist[0];
    for(int i=1;i<sizea;i++){
      averageRatio=(averageRatio+ratiolist[i])/2;
    }

    return averageRatio;
  }

  int Partitioner::getMacroStdcellAreaRatio(std::shared_ptr<PlacerBase> &pb_){
    int stdcellAverageArea=0;
    int macroAverageArea=0;
    int technum = pb_->techs().size();
    std::vector<Technology*> techs=pb_->techs();
    for(auto tech : techs){
      int libcellnum=tech->libCells().size();
      for(int i=0;i<libcellnum;i++){
        LibCell* libcell=tech->libCells()[i];
        if(libcell->isMacro()){
          macroAverageArea=(macroAverageArea+libcell->sizeX()*libcell->sizeY())/2;
        }else{
          stdcellAverageArea=(stdcellAverageArea+libcell->sizeX()*libcell->sizeY())/2;
        }
      }
    }
    int ratio=macroAverageArea/stdcellAverageArea;
    return ratio;
  }
}
