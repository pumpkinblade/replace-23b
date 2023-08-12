#include "parser.h"
#include "placerBase.h"
#include "log.h"
#include <fstream>
#include <memory>

namespace replace
{
  struct LIBPIN_DESC
  {
    std::string name;
    int x, y;
  };

  struct LIBCELL_DESC
  {
    std::string name;
    bool isMacro;
    int sizeX, sizeY;
    std::vector<LIBPIN_DESC> pins;
  };

  struct TECH_DESC
  {
    std::string name;
    std::vector<LIBCELL_DESC> cells;
  };

  // INST_DESC include instance name and used cell name
  struct INST_DESC
  {
    std::string name;
    std::string cellName;
  };

  struct PIN_DESC
  {
    std::string instName;
    std::string pinName;
  };

  struct NET_DESC
  {
    std::string name;
    std::vector<PIN_DESC> pins;
  };

  struct DIE_DESC
  {
    int lx, ly, ux, uy;
    int maxUtil;
    int rowStartX, rowStartY;
    int rowSizeX, rowSizeY;
    int repeatCount;
    std::string techName;
  };

  struct PLACER_23B_DESC
  {
    std::vector<TECH_DESC> techs;
    std::vector<INST_DESC> insts;
    std::vector<NET_DESC> nets;
    DIE_DESC topDieDesc;
    DIE_DESC bottomDieDesc;
    int termSizeX, termSizeY;
    int termSpace;
    int termCost;
  };

  void txtToDesc(const std::string &txtFilename, PLACER_23B_DESC *desc)
  {
    std::ifstream in(txtFilename);
    std::string token;
    if(in.bad())
    {
      LOG_CRITICAL("Couldn't open file `{}`", txtFilename);
      exit(-1);
    }

    while (!(in >> token).eof())
    {
      if (token == "NumTechnologies")
      {
        TECH_DESC *currTech;
        LIBCELL_DESC *currCell;
        LIBPIN_DESC *currLibPin;
        size_t numTechs;
        size_t numCells;
        size_t numPins;
        char isMacro;

        in >> numTechs;
        desc->techs.resize(numTechs);
        for (size_t i = 0; i < numTechs; i++)
        {
          currTech = &desc->techs[i];
          in >> token >> currTech->name >> numCells;
          currTech->cells.resize(numCells);
          for (size_t j = 0; j < numCells; j++)
          {
            currCell = &currTech->cells[j];
            in >> token >> isMacro >> currCell->name >> currCell->sizeX >> currCell->sizeY >> numPins;
            currCell->isMacro = isMacro == 'Y';
            currCell->pins.resize(numPins);
            for (size_t k = 0; k < numPins; k++)
            {
              currLibPin = &currCell->pins[k];
              in >> token >> currLibPin->name >> currLibPin->x >> currLibPin->y;
            }
          }
        }
      }
      else if (token == "DieSize")
      {
        int lx, ly, ux, uy;
        in >> lx >> ly >> ux >> uy;
        desc->bottomDieDesc.lx = desc->topDieDesc.lx = lx;
        desc->bottomDieDesc.ly = desc->topDieDesc.ly = ly;
        desc->bottomDieDesc.ux = desc->topDieDesc.ux = ux;
        desc->bottomDieDesc.uy = desc->topDieDesc.uy = uy;
      }
      else if (token == "TopDieMaxUtil")
      {
        in >> desc->topDieDesc.maxUtil;
      }
      else if (token == "BottomDieMaxUtil")
      {
        in >> desc->bottomDieDesc.maxUtil;
      }
      else if (token == "TopDieRows")
      {
        in >> desc->topDieDesc.rowStartX >> desc->topDieDesc.rowStartY >> desc->topDieDesc.rowSizeX >> desc->topDieDesc.rowSizeY >> desc->topDieDesc.repeatCount;
      }
      else if (token == "BottomDieRows")
      {
        in >> desc->bottomDieDesc.rowStartX >> desc->bottomDieDesc.rowStartY >> desc->bottomDieDesc.rowSizeX >> desc->bottomDieDesc.rowSizeY >> desc->bottomDieDesc.repeatCount;
      }
      else if (token == "TopDieTech")
      {
        in >> desc->topDieDesc.techName;
      }
      else if (token == "BottomDieTech")
      {
        in >> desc->bottomDieDesc.techName;
      }
      else if (token == "TerminalSize")
      {
        in >> desc->termSizeX >> desc->termSizeY;
      }
      else if (token == "TerminalSpacing")
      {
        in >> desc->termSpace;
      }
      else if (token == "TerminalCost")
      {
        in >> desc->termCost;
      }
      else if (token == "NumInstances")
      {
        size_t numInsts;
        INST_DESC *currInst;

        in >> numInsts;
        desc->insts.resize(numInsts);
        for (size_t i = 0; i < numInsts; i++)
        {
          currInst = &desc->insts[i];
          in >> token >> currInst->name >> currInst->cellName;
        }
      }
      else if (token == "NumNets")
      {
        size_t numNets;
        size_t numPins;
        size_t pos;
        NET_DESC *currNet;
        PIN_DESC *currPin;
        std::string pindesc;

        in >> numNets;
        desc->nets.resize(numNets);
        for (size_t i = 0; i < numNets; i++)
        {
          currNet = &desc->nets[i];
          in >> token >> currNet->name >> numPins;
          currNet->pins.resize(numPins);
          for (size_t j = 0; j < numPins; j++)
          {
            currPin = &currNet->pins[j];
            in >> token >> pindesc;
            pos = pindesc.find_first_of('/');
            currPin->instName.assign(pindesc.begin(), pindesc.begin() + pos);
            currPin->pinName.assign(pindesc.begin() + pos + 1, pindesc.end());
          }
        }
      }
    }
  }

  void printDebugInfo(const PLACER_23B_DESC *desc)
  {
    size_t numTechs = desc->techs.size();
    size_t numCells;
    size_t numLibPins;
    size_t numInsts;
    size_t numNets;
    size_t numPins;
    const TECH_DESC *currTech;
    const LIBCELL_DESC *currCell;
    const LIBPIN_DESC *currLibPin;
    const INST_DESC *currInst;
    const PIN_DESC *currPin;
    const NET_DESC *currNet;

    LOG_DEBUG("NumTechnologies {}", numTechs);
    for (size_t i = 0; i < numTechs; i++)
    {
      currTech = &desc->techs[i];
      numCells = currTech->cells.size();
      LOG_DEBUG("Tech {} {}", currTech->name, numCells);
      for (size_t j = 0; j < numCells; j++)
      {
        currCell = &currTech->cells[j];
        numLibPins = currCell->pins.size();
        LOG_DEBUG("LibCell {} {} {} {} {}", currCell->isMacro, currCell->name,
                 currCell->sizeX, currCell->sizeY, numLibPins);
        for (size_t k = 0; k < numLibPins; k++)
        {
          currLibPin = &currCell->pins[k];
          LOG_DEBUG("Pin {} {} {}", currLibPin->name, currLibPin->x, currLibPin->y);
        }
      }
    }

    LOG_DEBUG("DieSize {} {} {} {}", desc->topDieDesc.lx, desc->topDieDesc.ly,
             desc->topDieDesc.ux, desc->topDieDesc.uy);
    LOG_DEBUG("TopDieMaxUtil {}", desc->topDieDesc.maxUtil);
    LOG_DEBUG("BottomDieMaxUtil {}", desc->bottomDieDesc.maxUtil);
    LOG_DEBUG("TopDieRows {} {} {} {} {}", desc->topDieDesc.rowStartX, desc->topDieDesc.rowStartY,
             desc->topDieDesc.rowSizeX, desc->topDieDesc.rowSizeY, desc->topDieDesc.repeatCount);
    LOG_DEBUG("BottomDieRows {} {} {} {} {}", desc->bottomDieDesc.rowStartX, desc->bottomDieDesc.rowStartY,
             desc->bottomDieDesc.rowSizeX, desc->bottomDieDesc.rowSizeY, desc->bottomDieDesc.repeatCount);
    LOG_DEBUG("TopDieTech {}", desc->topDieDesc.techName);
    LOG_DEBUG("BottomDieTech {}", desc->bottomDieDesc.techName);
    LOG_DEBUG("TerminalSize {} {}", desc->termSizeX, desc->termSizeY);
    LOG_DEBUG("TerminalSpacing {}", desc->termSpace);
    LOG_DEBUG("TerminalCost {}", desc->termCost);

    numInsts = desc->insts.size();
    LOG_DEBUG("NumInstances {}", numInsts);
    for (size_t i = 0; i < numInsts; i++)
    {
      currInst = &desc->insts[i];
      LOG_DEBUG("Inst {} {}", currInst->name, currInst->cellName);
    }

    numNets = desc->nets.size();
    LOG_DEBUG("NumNets {}", numNets);
    for (size_t i = 0; i < numNets; i++)
    {
      currNet = &desc->nets[i];
      numPins = currNet->pins.size();
      LOG_DEBUG("Net {} {}", currNet->name, numPins);
      for (size_t j = 0; j < numPins; j++)
      {
        currPin = &currNet->pins[j];
        LOG_DEBUG("Pin {}/{}", currPin->instName, currPin->pinName);
      }
    }
  }

  std::shared_ptr<PlacerBase> Parser::txtToPlacerBase(const std::string &txtFilename)
  {
    // auto pb = std::make_shared<Placer23b>();
    auto pb = std::make_shared<PlacerBase>();

    PLACER_23B_DESC desc;
    LOG_TRACE("Parse txt file begin");
    txtToDesc(txtFilename, &desc);
    LOG_TRACE("Parse txt file end");

    // Init libCellId & libPinId
    std::unordered_map<std::string, int> libCellNameIdMap;
    std::vector<std::unordered_map<std::string, int>> libPinNameIdMaps;
    {
      TECH_DESC* techDesc = &desc.techs.front();
      libCellNameIdMap.reserve(techDesc->cells.size());
      libPinNameIdMaps.reserve(techDesc->cells.size());
      for(int i = 0; i < techDesc->cells.size(); i++)
      {
        LIBCELL_DESC* cellDesc = &techDesc->cells[i];
        libCellNameIdMap.emplace(cellDesc->name, i);
        libPinNameIdMaps.emplace_back();
        for(int j = 0; j < cellDesc->pins.size(); j++)
        {
          libPinNameIdMaps[i].emplace(cellDesc->pins[j].name, j);
        }
      }
    }

    // Process technology
    pb->techStor_.reserve(desc.techs.size());
    for (size_t i = 0; i < desc.techs.size(); i++)
    {
      pb->techStor_.emplace_back(desc.techs[i].name);
      Technology* tech = &pb->techStor_.back();
      pb->techs_.push_back(tech);
      pb->techNameMap_.emplace(tech->name(), tech);

      TECH_DESC *techDesc = &desc.techs[i];
      size_t numPins = 0;
      size_t numCells = techDesc->cells.size();
      for (const auto &cell : desc.techs[i].cells)
        numPins += cell.pins.size();
      tech->cellStor_.reserve(numCells);
      tech->cells_.resize(numCells);
      tech->pinStor_.reserve(numPins);

      for (size_t j = 0; j < numCells; j++)
      {
        LIBCELL_DESC *cellDesc = &techDesc->cells[j];
        int cellId = libCellNameIdMap.at(cellDesc->name);
        tech->cellStor_.emplace_back(cellId, cellDesc->sizeX, cellDesc->sizeY, cellDesc->isMacro);
        LibCell *cell = &tech->cellStor_.back();
        cell->pins_.resize(cellDesc->pins.size());

        tech->cells_[cellId] = cell;

        for (const auto &pinDesc : cellDesc->pins)
        {
          int pinId = libPinNameIdMaps[cellId].at(pinDesc.name);
          tech->pinStor_.emplace_back(pinId, pinDesc.x, pinDesc.y);
          cell->pins_[pinId] = &tech->pinStor_.back();
        }
      }
    }

    // Process dies
    pb->dieStor_.reserve(3);
    // top die
    pb->dieStor_.emplace_back();
    pb->dieStor_.back().setName("top");
    pb->dieStor_.back().setDieBox(desc.topDieDesc.lx, desc.topDieDesc.ly,
                                  desc.topDieDesc.ux, desc.topDieDesc.uy);
    pb->dieStor_.back().setRowParams(desc.topDieDesc.rowStartX, desc.topDieDesc.rowStartY,
                                     desc.topDieDesc.rowSizeX, desc.topDieDesc.rowSizeY,
                                     desc.topDieDesc.repeatCount);
    pb->dieStor_.back().setMaxUtil(desc.topDieDesc.maxUtil / 100.f);
    pb->dieStor_.back().setTech(pb->tech(desc.topDieDesc.techName));
    pb->dieNameMap_.emplace("top", &pb->dieStor_.back());
    // bottom die
    pb->dieStor_.emplace_back();
    pb->dieStor_.back().setName("bottom");
    pb->dieStor_.back().setDieBox(desc.bottomDieDesc.lx, desc.bottomDieDesc.ly,
                                  desc.bottomDieDesc.ux, desc.bottomDieDesc.uy);
    pb->dieStor_.back().setRowParams(desc.bottomDieDesc.rowStartX, desc.bottomDieDesc.rowStartY,
                                     desc.bottomDieDesc.rowSizeX, desc.bottomDieDesc.rowSizeY,
                                     desc.bottomDieDesc.repeatCount);
    pb->dieStor_.back().setMaxUtil(desc.bottomDieDesc.maxUtil / 100.f);
    pb->dieStor_.back().setTech(pb->tech(desc.bottomDieDesc.techName));
    pb->dieNameMap_.emplace("bottom", &pb->dieStor_.back());
    // terminal die
    pb->dieStor_.emplace_back();
    pb->dieStor_.back().setName("terminal");
    pb->dieStor_.back().setDieBox(desc.topDieDesc.lx, desc.topDieDesc.ly,
                                  desc.topDieDesc.ux, desc.topDieDesc.uy);
    pb->dieStor_.back().setCoreBox(desc.topDieDesc.lx, desc.topDieDesc.ly,
                                   desc.topDieDesc.ux, desc.topDieDesc.uy);
    pb->dieNameMap_.emplace("terminal", &pb->dieStor_.back());
    // Note that if a die is in pb->dies_, then it will participate in global placement
    // therefore, we only push top die into pb->dies_
    pb->dies_.push_back(pb->die("top"));
    pb->dies_.push_back(pb->die("bottom"));
    pb->dies_.push_back(pb->die("terminal"));

    // Process terminal
    pb->termSizeX_ = desc.termSizeX;
    pb->termSizeY_ = desc.termSizeY;
    pb->termSpace_ = desc.termSpace;
    pb->termCost_ = desc.termCost;

    // Process netlist
    pb->instStor_.reserve(desc.insts.size());
    pb->insts_.reserve(desc.insts.size());
    pb->netStor_.reserve(desc.nets.size());
    pb->nets_.reserve(desc.nets.size());
    pb->netNameStor_.reserve(desc.nets.size());
    size_t numPins = 0;
    for(const auto& netDesc : desc.nets)
    {
      numPins += netDesc.pins.size();
    }
    pb->pinStor_.reserve(numPins);
    pb->pins_.reserve(numPins);

    // Process inst
    std::unordered_map<std::string, Instance*> instNameMap;
    for (const auto& instDesc : desc.insts)
    {
      pb->instStor_.emplace_back();
      Instance* inst = &pb->instStor_.back();
      pb->insts_.push_back(inst);
      pb->die("top")->addInstance(inst);
      instNameMap.emplace(instDesc.name, inst);

      // find libcell
      // assume than all insts are on top die
      int libcellId = libCellNameIdMap[instDesc.cellName];
      LibCell* libcell = pb->die("top")->tech()->libCell(libcellId);

      // inst init
      inst->setName(instDesc.name);
      inst->setLibCellId(libcellId);
      inst->setFixed(false);
      inst->setMacro(libcell->isMacro());
      int lx = pb->die("top")->coreLx();
      int ly = pb->die("top")->coreLy();
      inst->setBox(lx, ly, lx + libcell->sizeX(), ly + libcell->sizeY());
    }
    
    // Process net
    for (const auto& netDesc : desc.nets)
    {
      pb->netStor_.emplace_back();
      pb->netNameStor_.emplace_back(netDesc.name);
      Net* net = &pb->netStor_.back();
      pb->nets_.push_back(net);

      for (const auto& pinDesc : netDesc.pins)
      {
        // find instance
        Instance* inst = instNameMap.at(pinDesc.instName);
        // find libcell & libpin
        int libCellId = inst->libCellId();
        int libPinId = libPinNameIdMaps[libCellId].at(pinDesc.pinName);
        LibCell* libcell = pb->die("top")->tech()->libCell(libCellId);
        LibPin* libpin = libcell->libPin(libPinId);

        // add pin
        pb->pinStor_.emplace_back();
        Pin* pin = &pb->pinStor_.back();
        pb->pins_.push_back(pin);
        pin->setLibPinId(libPinId);
        pin->setInstance(inst);
        pin->setNet(net);
        pin->updateLocation(inst, libpin->x(), libpin->y());

        // assign inst & net
        net->addPin(pin);
        inst->addPin(pin);
      }
      net->updateBox();
    }

    return pb;
  }
}
