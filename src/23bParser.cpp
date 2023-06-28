#include "parser.h"
#include "log.h"
#include <fstream>

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
    size_t repeatCount;
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

  void ParseTxt(const std::string& txtFilename, PLACER_23B_DESC* desc)
  {
    std::ifstream in(txtFilename);
    std::string token;

    while(!(in >> token).eof())
    {
      if(token == "NumTechnologies")
      {
        TECH_DESC* currTech;
        LIBCELL_DESC* currCell;
        LIBPIN_DESC* currLibPin;
        size_t numTechs;
        size_t numCells;
        size_t numPins;
        char isMacro;

        in >> numTechs;
        desc->techs.resize(numTechs);
        for(size_t i = 0; i < numTechs; i++)
        {
          currTech = &desc->techs[i];
          in >> token >> currTech->name >> numCells;
          currTech->cells.resize(numCells);
          for(size_t j = 0; j < numCells; j++)
          {
            currCell = &currTech->cells[j];
            in >> token >> isMacro >> currCell->name 
               >> currCell->sizeX >> currCell->sizeY >> numPins;
            currCell->isMacro = isMacro == 'Y';
            currCell->pins.resize(numPins);
            for(size_t k = 0; k < numPins; k++)
            {
              currLibPin = &currCell->pins[k];
              in >> token >> currLibPin->name >> currLibPin->x >> currLibPin->y;
            }
          }
        }
      }
      else if(token == "DieSize")
      {
        int lx, ly, ux, uy;
        in >> lx >> ly >> ux >> uy;
        desc->bottomDieDesc.lx = desc->topDieDesc.lx = lx;
        desc->bottomDieDesc.ly = desc->topDieDesc.ly = ly;
        desc->bottomDieDesc.ux = desc->topDieDesc.ux = ux;
        desc->bottomDieDesc.uy = desc->topDieDesc.uy = uy;
      }
      else if(token == "TopDieMaxUtil")
      {
        in >> desc->topDieDesc.maxUtil;
      }
      else if(token == "BottomDieMaxUtil")
      {
        in >> desc->bottomDieDesc.maxUtil;
      }
      else if(token == "TopDieRows")
      {
        in >> desc->topDieDesc.rowStartX >> desc->topDieDesc.rowStartY
           >> desc->topDieDesc.rowSizeX >> desc->topDieDesc.rowSizeY
           >> desc->topDieDesc.repeatCount;
      }
      else if(token == "BottomDieRows")
      {
        in >> desc->bottomDieDesc.rowStartX >> desc->bottomDieDesc.rowStartY
           >> desc->bottomDieDesc.rowSizeX >> desc->bottomDieDesc.rowSizeY
           >> desc->bottomDieDesc.repeatCount;
      }
      else if(token == "TopDieTech")
      {
        in >> desc->topDieDesc.techName;
      }
      else if(token == "BottomDieTech")
      {
        in >> desc->bottomDieDesc.techName;
      }
      else if(token == "TerminalSize")
      {
        in >> desc->termSizeX >> desc->termSizeY;
      }
      else if(token == "TerminalSpacing")
      {
        in >> desc->termSpace;
      }
      else if(token == "TerminalCost")
      {
        in >> desc->termCost;
      }
      else if(token == "NumInstances")
      {
        size_t numInsts;
        INST_DESC* currInst;

        in >> numInsts;
        desc->insts.resize(numInsts);
        for(size_t i = 0; i < numInsts; i++)
        {
          currInst = &desc->insts[i];
          in >> token >> currInst->name >> currInst->cellName;
        }
      }
      else if(token == "NumNets")
      {
        size_t numNets;
        size_t numPins;
        size_t pos;
        NET_DESC* currNet;
        PIN_DESC* currPin;
        std::string pindesc;

        in >> numNets;
        desc->nets.resize(numNets);
        for(size_t i = 0; i < numNets; i++)
        {
          currNet = &desc->nets[i];
          in >> token >> currNet->name >> numPins;
          currNet->pins.resize(numPins);
          for(size_t j = 0; j < numPins; j++)
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

  void PrintInfo(PLACER_23B_DESC* desc)
  {
    size_t numTechs = desc->techs.size();
    size_t numCells;
    size_t numLibPins;
    size_t numInsts;
    size_t numNets;
    size_t numPins;
    TECH_DESC* currTech;
    LIBCELL_DESC* currCell;
    LIBPIN_DESC* currLibPin;
    INST_DESC* currInst;
    PIN_DESC* currPin;
    NET_DESC* currNet;

    LOG_INFO("NumTechnologies {}", numTechs);
    for(size_t i = 0; i < numTechs; i++)
    {
      currTech = &desc->techs[i];
      numCells = currTech->cells.size();
      LOG_INFO("Tech {} {}", currTech->name, numCells);
      for(size_t j = 0; j < numCells; j++)
      {
        currCell = &currTech->cells[j];
        numLibPins = currCell->pins.size();
        LOG_INFO("LibCell {} {} {} {} {}", currCell->isMacro, currCell->name,
                 currCell->sizeX, currCell->sizeY, numLibPins);
        for(size_t k = 0; k < numLibPins; k++)
        {
          currLibPin = &currCell->pins[k];
          LOG_INFO("Pin {} {} {}", currLibPin->name, currLibPin->x, currLibPin->y);
        }
      }
    }

    LOG_INFO("DieSize {} {} {} {}", desc->topDieDesc.lx, desc->topDieDesc.ly,
             desc->topDieDesc.ux, desc->topDieDesc.uy);
    LOG_INFO("TopDieMaxUtil {}", desc->topDieDesc.maxUtil);
    LOG_INFO("BottomDieMaxUtil {}", desc->bottomDieDesc.maxUtil);
    LOG_INFO("TopDieRows {} {} {} {} {}", desc->topDieDesc.rowStartX, desc->topDieDesc.rowStartY,
             desc->topDieDesc.rowSizeX, desc->topDieDesc.rowSizeY, desc->topDieDesc.repeatCount);
    LOG_INFO("BottomDieRows {} {} {} {} {}", desc->bottomDieDesc.rowStartX, desc->bottomDieDesc.rowStartY,
             desc->bottomDieDesc.rowSizeX, desc->bottomDieDesc.rowSizeY, desc->bottomDieDesc.repeatCount);
    LOG_INFO("TopDieTech {}", desc->topDieDesc.techName);
    LOG_INFO("BottomDieTech {}", desc->bottomDieDesc.techName);
    LOG_INFO("TerminalSize {} {}", desc->termSizeX, desc->termSizeY);
    LOG_INFO("TerminalSpacing {}", desc->termSpace);
    LOG_INFO("TerminalCost {}", desc->termCost);

    numInsts = desc->insts.size();
    LOG_INFO("NumInstances {}", numInsts);
    for(size_t i = 0; i < numInsts; i++)
    {
      currInst = &desc->insts[i];
      LOG_INFO("Inst {} {}", currInst->name, currInst->cellName);
    }

    numNets = desc->nets.size();
    LOG_INFO("NumNets {}", numNets);
    for(size_t i = 0; i < numNets; i++)
    {
      currNet = &desc->nets[i];
      numPins = currNet->pins.size();
      LOG_INFO("Net {} {}", currNet->name, numPins);
      for(size_t j = 0; j < numPins; j++)
      {
        currPin = &currNet->pins[j];
        LOG_INFO("Pin {}/{}", currPin->instName, currPin->pinName);
      }
    }
  }

  std::shared_ptr<Placer23b> Parser::TxtToPlacer23b(const std::string& txtFilename)
  {
    auto pb = std::make_shared<Placer23b>();

    PLACER_23B_DESC desc;
    LOG_TRACE("Parse txt file begin");
    ParseTxt(txtFilename, &desc);
    LOG_TRACE("Parse txt file end");

    // Process technology
    pb->techStor_.reserve(desc.techs.size());
    for(size_t i = 0; i < desc.techs.size(); i++)
    {
      pb->techStor_.emplace_back();
      Technology* tech = &pb->techStor_.back();
      TECH_DESC* techDesc = &desc.techs[i];
      size_t numPins = 0;
      size_t numCells = techDesc->cells.size();
      for(const auto& cell : desc.techs[i].cells)
        numPins += cell.pins.size();
      pb->techStor_[i].cellStor_.reserve(numCells);
      pb->techStor_[i].pinStor_.reserve(numPins);
      pb->techNameMap_.emplace(techDesc->name, tech);

      for(size_t j = 0; j < numCells; j++)
      {
        LIBCELL_DESC* cellDesc = &techDesc->cells[j];
        tech->cellStor_.emplace_back(cellDesc->sizeX, cellDesc->sizeY, cellDesc->isMacro);
        LibCell* cell = &tech->cellStor_.back();

        tech->cellNameMap_.emplace(cellDesc->name, cell);
        tech->cells_.push_back(cell);
        if(cell->isMacro())
          tech->macros_.push_back(cell);
        else
          tech->stdCells_.push_back(cell);
        tech->siteSizeX_ = tech->siteSizeY_ = 0;

        for(size_t k = 0; k < cellDesc->pins.size(); k++)
        {
          tech->pinStor_.emplace_back(cellDesc->pins[k].x, cellDesc->pins[k].y);
          cell->addPin(cellDesc->pins[k].name, &tech->pinStor_.back());
        }
      }
    }

    // Process die
    pb->topDie_.setDieBox(desc.topDieDesc.lx, desc.topDieDesc.ly,
                          desc.topDieDesc.ux, desc.topDieDesc.uy);
    pb->bottomDie_.setDieBox(desc.bottomDieDesc.lx, desc.bottomDieDesc.ly,
                             desc.bottomDieDesc.ux, desc.bottomDieDesc.uy);
    for(size_t i = 0; i < desc.topDieDesc.repeatCount; i++)
    {
      pb->topDie_.addRow(Row(desc.topDieDesc.rowStartX, 
                             desc.topDieDesc.rowStartY + i*desc.topDieDesc.rowSizeY,
                             desc.topDieDesc.rowStartX + desc.topDieDesc.rowSizeX,
                             desc.topDieDesc.rowStartY + (i+1)*desc.topDieDesc.rowSizeY));
    }
    pb->topDie_.updateCoreBox();
    for(size_t i = 0; i < desc.bottomDieDesc.repeatCount; i++)
    {
      pb->topDie_.addRow(Row(desc.bottomDieDesc.rowStartX, 
                             desc.bottomDieDesc.rowStartY + i*desc.bottomDieDesc.rowSizeY,
                             desc.bottomDieDesc.rowStartX + desc.bottomDieDesc.rowSizeX,
                             desc.bottomDieDesc.rowStartY + (i+1)*desc.bottomDieDesc.rowSizeY));
    }
    pb->bottomDie_.updateCoreBox();
    pb->topUtil_ = desc.topDieDesc.maxUtil / 100.f;
    pb->bottomUtil_ = desc.bottomDieDesc.maxUtil / 100.f;
    pb->topTech_ = pb->techNameMap_.at(desc.topDieDesc.techName);
    pb->bottomTech_ = pb->techNameMap_.at(desc.bottomDieDesc.techName);

    // Process term
    pb->termSizeX_ = desc.termSizeX;
    pb->termSizeY_ = desc.termSizeY;
    pb->termSpace_ = desc.termSpace;
    pb->termCost_ = desc.termCost;

    // Process netlist
    pb->instStor_.reserve(desc.insts.size());
    pb->netStor_.reserve(desc.nets.size());
    size_t numPins = 0;
    for(const auto& n : desc.nets)
      numPins += n.pins.size();
    pb->pinStor_.reserve(numPins);
    // Process inst
    for(size_t i = 0; i < desc.insts.size(); i++)
    {
      pb->instStor_.emplace_back(desc.insts[i].cellName);
      pb->instNameMap_.emplace(desc.insts[i].name, &pb->instStor_.back());
      pb->insts_.push_back(&pb->instStor_.back());
    }
    // Process net and pin
    for(size_t i = 0; i < desc.nets.size(); i++)
    {
      pb->netStor_.emplace_back();
      pb->netNameMap_.emplace(desc.nets[i].name, &pb->netStor_.back());
      pb->nets_.push_back(&pb->netStor_.back());
      for(size_t j = 0; j < desc.nets[i].pins.size(); j++)
      {
        pb->pinStor_.emplace_back(desc.nets[i].pins[j].instName,
                                  desc.nets[i].pins[j].pinName);
        pb->pins_.push_back(&pb->pinStor_.back());
        pb->netStor_.back().addPin(&pb->pinStor_.back());
      }
    }

    return pb;
  }
}