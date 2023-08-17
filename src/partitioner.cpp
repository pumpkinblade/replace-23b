#include "technology.h"
#include "log.h"
#include "placerBase.h"
#include "partitioner.h"
#include "replace.h"
#include <unordered_map>
#include <algorithm>
#include <unordered_set>
#include <libkahypar.h>
#include <thread>
#include <iostream>

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
    LOG_DEBUG("partition macros, top: {} bottom: {}",
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

  // void Partitioner::hmetistest(std::shared_ptr<PlacerBase> pb_){
  //   // 将pb中的instance和net转换为hmetis的数据结构
  //   // 1. instance转换为hmetis的hypervertex
  //   // 2. net转换为hmetis的edge
  //   // 3. 调用hmetis进行partition
  //   LOG_INFO("start hmetis partition");
  //   // 获得instance数量
  //   int nvtxs = pb_->insts().size();
  //   // 获得net数量
  //   int nhedges = pb_->nets().size();
  //   // 获得instance的weight, 目前暂时设为1
  //   int *vwgts = new int[nvtxs];
  //   // 将net转换为eptr和eind
  //   int *eptr = new int[nhedges + 1];
  //   int eindnum = 0;
  //   for(auto net : pb_->nets()){
  //     eindnum += net->pins().size();
  //   }
  //   int *eind = new int[eindnum];
  //   int i=0;
  //   for(auto net : pb_->nets()){
  //     for(auto pin : net->pins()){
  //       eind[i] = pin->instance()->extId();
  //       i++;
  //     }
  //   }
  //   // 获得hmetis的参数
  //   int *options = new int[10];
  //   options[0] = 0;
  //   int *part = new int[nvtxs];
  //   int *edgecut = new int[1];
  //   int nparts = 2;
  //   int ubfactor = 10;
  //   int *hewgts = new int[nhedges];
  //   for(int i=0; i<nhedges; i++){
  //     hewgts[i] = 1;
  //   }
  //   // 调用hmetis进行partition
  //   HMETIS_PartRecursive(nvtxs, nhedges, NULL, *eptr, *eind, NULL, nparts, ubfactor, *options, *part, *edgecut);
  //   std::string filename = "hmetis.txt";
  //   std::ofstream out(filename);
  //   for(int i=0; i<nvtxs; i++){
  //     out << i << " " << part[i] << std::endl;
  //   }
  //   LOG_INFO("hmetis partition finished");
  // }

  void moveDecide(int64_t topArea, int64_t botArea, int64_t topCap, int64_t botCap,
                  bool* moveToTop, bool* moveToBot)
  {
    bool topCanContain = topArea < topCap;
    bool botCanContain = botArea < botCap;

    *moveToTop = *moveToBot = false;
    if(topCanContain && botCanContain)
    {
      *moveToBot = ((float)std::rand() / RAND_MAX) < 0.5;
      *moveToTop = !(*moveToBot);
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
        LibCell* libcell = botdie->tech()->libCell(inst->libCellName());
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
        LibCell* libcell = botdie->tech()->libCell(inst->libCellName());
        int instTopArea, instBotArea;
        instTopArea = inst->dx() * inst->dy();
        instBotArea = libcell->sizeX() * libcell->sizeY();

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

  char const *const ini_file =
          "/home/workspace/replace-23bmain/replace-23b/test/cut_rKaHyPar_sea20.ini";


  // Adapted example from the README to have something to hammer on
  auto run_kahypar() {
    kahypar_context_t *context = kahypar_context_new();
    kahypar_configure_context_from_file(context, ini_file);

    const kahypar_hypernode_id_t num_vertices = 7;
    const kahypar_hyperedge_id_t num_hyperedges = 4;

    std::unique_ptr<kahypar_hyperedge_weight_t[]> hyperedge_weights =
            std::make_unique<kahypar_hyperedge_weight_t[]>(4);

    // force the cut to contain hyperedge 0 and 2
    hyperedge_weights[0] = 1;
    hyperedge_weights[1] = 1000;
    hyperedge_weights[2] = 1;
    hyperedge_weights[3] = 1000;

    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);

    hyperedge_indices[0] = 0;
    hyperedge_indices[1] = 2;
    hyperedge_indices[2] = 6;
    hyperedge_indices[3] = 9;
    hyperedge_indices[4] = 12;

    std::unique_ptr<kahypar_hyperedge_id_t[]> hyperedges =
            std::make_unique<kahypar_hyperedge_id_t[]>(12);

    // hypergraph from hMetis manual page 14
    hyperedges[0] = 0;
    hyperedges[1] = 2;
    hyperedges[2] = 0;
    hyperedges[3] = 1;
    hyperedges[4] = 3;
    hyperedges[5] = 4;
    hyperedges[6] = 3;
    hyperedges[7] = 4;
    hyperedges[8] = 6;
    hyperedges[9] = 2;
    hyperedges[10] = 5;
    hyperedges[11] = 6;

    const double imbalance = 0.03;
    const kahypar_partition_id_t k = 2;

    kahypar_hyperedge_weight_t objective = 0;

    std::vector<kahypar_partition_id_t> partition(num_vertices, -1);

    kahypar_partition(num_vertices, num_hyperedges, imbalance, k,
                      /*vertex_weights */ nullptr, hyperedge_weights.get(),
                      hyperedge_indices.get(), hyperedges.get(), &objective,
                      context, partition.data());
    for(int i = 0; i != num_vertices; ++i) {
      std::cout << i << ":" << partition[i] << std::endl;
    }

    kahypar_context_free(context);
  }

  void Partitioner::do_run_kahypar(){
    // auto const task = [&] {
    //     while (true) { run_kahypar(); }
    // };

    // auto threads = std::vector<std::thread>{};
    // for (int i = 0; i < std::thread::hardware_concurrency(); ++i) {
    //     threads.emplace_back(task);
    // }
    // for (auto &thread: threads) { thread.join(); }
    kahypar_context_t* context = kahypar_context_new();
    kahypar_configure_context_from_file(context, "/home/workspace/replace-23bmain/replace-23b/test/cut_rKaHyPar_sea20.ini");
    
    kahypar_set_seed(context, 42);

    const kahypar_hypernode_id_t num_vertices = 7;
    const kahypar_hyperedge_id_t num_hyperedges = 4;

    std::unique_ptr<kahypar_hyperedge_weight_t[]> hyperedge_weights = std::make_unique<kahypar_hyperedge_weight_t[]>(4);

    // force the cut to contain hyperedge 0 and 2
    hyperedge_weights[0] = 1;  hyperedge_weights[1] = 1000;
    hyperedge_weights[2] = 1;  hyperedge_weights[3] = 1000;

    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(5);

    hyperedge_indices[0] = 0; hyperedge_indices[1] = 2;
    hyperedge_indices[2] = 6; hyperedge_indices[3] = 9;
    hyperedge_indices[4] = 12;

    std::unique_ptr<kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<kahypar_hyperedge_id_t[]>(12);

    // hypergraph from hMetis manual page 14
    hyperedges[0] = 0;  hyperedges[1] = 2;
    hyperedges[2] = 0;  hyperedges[3] = 1;
    hyperedges[4] = 3;  hyperedges[5] = 4;
    hyperedges[6] = 3;  hyperedges[7] = 4;
    hyperedges[8] = 6;  hyperedges[9] = 2;
    hyperedges[10] = 5; hyperedges[11] = 6;

    const double imbalance = 0.01;
    const kahypar_partition_id_t k = 2;

    kahypar_hyperedge_weight_t objective = 0;

    std::vector<kahypar_partition_id_t> partition(num_vertices, -1);

    kahypar_partition(num_vertices, num_hyperedges,
                      imbalance, k,
                      /*vertex_weights */ nullptr, hyperedge_weights.get(),
                      hyperedge_indices.get(), hyperedges.get(),
                      &objective, context, partition.data());

    for(int i = 0; i != num_vertices; ++i) {
      std::cout << i << ":" << partition[i] << std::endl;
    }

    kahypar_context_free(context);
  }

  void Partitioner::partitionInstance(std::shared_ptr<PlacerBase> &pb_){
    // convert the instance and net to hypergraph
    LOG_TRACE("start partition");
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
    kahypar_configure_context_from_file(context, "/home/workspace/replace-23bmain/replace-23b/test/cut_rKaHyPar_sea20.ini");
    
    kahypar_set_seed(context, 42);

    kahypar_partition(num_vertices, num_hyperedges,
                      imbalance, k,
                      /*vertex_weights */ nullptr, hyperedge_weights.get(),
                      hyperedge_indices.get(), hyperedges.get(),
                      &objective, context, partition.data());

    // for(int i = 0; i != num_vertices; ++i) {
    //   std::cout << i << ":" << partition[i] << std::endl;
    // }

    // move the instance to the bottom die
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
    for (int i=0; i<instanceNum; i++)
    {
      Instance *instance = pb_->insts()[i];
      // colloect all macros
      if(instance->isMacro()){
        macros.push_back(instance);
        continue;
      }
      bool moveToBottom = partition[i] ==1 ? true : false;
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
    kahypar_context_free(context);
    LOG_TRACE("finish partition");
  }


}
