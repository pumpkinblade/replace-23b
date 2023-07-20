#include "parser.h"
#include "log.h"
#include "technology.h"
#include <algorithm>
#include <lefrReader.hpp>
#include <defrReader.hpp>

namespace replace
{
  struct LEF_DATABASE
  {
    double lefUnit;
    lefiSite lefCoreSite;
    std::unordered_map<std::string, lefiMacro> lefMacroMap;
    std::unordered_map<std::string, int> lefMacroPinIdxMap;
    std::vector<std::vector<lefiPin>> lefMacroPins;
  };

  static int lefUnitCallback(lefrCallbackType_e typ, lefiUnits *unit, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefUnit = unit->databaseNumber();
    return 0;
  }

  static int lefSiteCallback(lefrCallbackType_e typ, lefiSite *site, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    if(std::strcmp(site->name(), "CoreSite") == 0)
      db->lefCoreSite = *site;
    return 0;
  }

  static int lefMacroCallback(lefrCallbackType_e typ, lefiMacro *macro, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefMacroMap.emplace(macro->name(), *macro);
    return 0;
  }

  static int lefMacroBeginCallback(lefrCallbackType_e typ, const char *name, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefMacroPinIdxMap.emplace(name, db->lefMacroPins.size());
    db->lefMacroPins.emplace_back();
    return 0;
  }

  static int lefPinCallback(lefrCallbackType_e typ, lefiPin *pin, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefMacroPins.back().push_back(*pin);
    return 0;
  }

  void parseLef(const std::string &lefFilename, LEF_DATABASE *lefdb)
  {
    // init lef reader and set callback
    lefrInit();
    lefrSetUnitsCbk(lefUnitCallback);
    lefrSetSiteCbk(lefSiteCallback);
    lefrSetMacroBeginCbk(lefMacroBeginCallback);
    lefrSetMacroCbk(lefMacroCallback);
    lefrSetPinCbk(lefPinCallback);

    // open the lef file
    FILE *lefFile = nullptr;
    if ((lefFile = fopen(lefFilename.c_str(), "r")) == 0)
    {
      LOG_CRITICAL("Couldn't open LEF File `{}`", lefFilename);
      exit(1);
    }

    // invoke the lef reader
    int res = lefrRead(lefFile, lefFilename.c_str(), lefdb);
    fclose(lefFile);
    if (res != 0)
    {
      LOG_CRITICAL("An error occurs when parsing LEF file `{}`", lefFilename);
      exit(1);
    }
  }

  std::pair<double, double> getPinLocation(const lefiPin &pin)
  {
    double lx = INT_MAX;
    double ux = INT_MIN;
    double ly = INT_MAX;
    double uy = INT_MIN;
    for (int i = 0; i < pin.numPorts(); i++)
    {
      lefiGeometries *geo = pin.port(i);
      for (int j = 0; j < geo->numItems(); j++)
      {
        if (geo->itemType(j) == lefiGeomRectE)
        {
          lx = std::min(lx, geo->getRect(j)->xl);
          ux = std::max(ux, geo->getRect(j)->xh);
          ly = std::min(ly, geo->getRect(j)->yl);
          uy = std::max(uy, geo->getRect(j)->yh);
        }
      }
    }
    return std::make_pair((lx + ux) / 2, (ly + uy) / 2);
  }

  std::shared_ptr<Technology> Parser::lefToTechnology(const std::string& lefFilename)
  {
    auto tech = std::make_shared<Technology>("lef");

    LEF_DATABASE lefdb;
    LOG_TRACE("Parse LEF Begin");
    parseLef(lefFilename, &lefdb);
    LOG_TRACE("Parse LEF End");

    LOG_TRACE("Process lef site");
    int siteX = static_cast<int>(lefdb.lefCoreSite.sizeX() * lefdb.lefUnit);
    int siteY = static_cast<int>(lefdb.lefCoreSite.sizeY() * lefdb.lefUnit);
    tech->siteSizeX_ = siteX;
    tech->siteSizeY_ = siteY;

    LOG_TRACE("Process lef libcell and libpin");
    tech->cellStor_.reserve(lefdb.lefMacroMap.size());
    size_t numPins = 0;
    for(const auto& p : lefdb.lefMacroPins)
      numPins += p.size();
    tech->pinStor_.reserve(numPins);

    // add libcell and libpin to tech
    for (const auto& it : lefdb.lefMacroMap)
    {
      int w = static_cast<int>(it.second.sizeX() * lefdb.lefUnit);
      int h = static_cast<int>(it.second.sizeY() * lefdb.lefUnit);
      bool isMacro = h > siteY;
      tech->cellStor_.emplace_back(it.first, w, h, isMacro);
      tech->cellNameMap_.emplace(it.first, &tech->cellStor_.back());
      tech->cells_.push_back(&tech->cellStor_.back());

      int idx = lefdb.lefMacroPinIdxMap[it.first];
      const auto& pins = lefdb.lefMacroPins[idx];
      for (const auto &pin : pins)
      {
        std::pair<double, double> offset = getPinLocation(pin);
        int x = static_cast<int>(offset.first * lefdb.lefUnit);
        int y = static_cast<int>(offset.second * lefdb.lefUnit);
        tech->pinStor_.emplace_back(pin.name(), x, y);
        tech->cellStor_.back().addLibPin(&tech->pinStor_.back());
      }
    }
    return tech;
  }

  struct DEF_DATABASE
  {
    defiBox defDieArea;
    double defUnit;
    std::vector<defiRow> defRows;
    std::unordered_map<std::string, defiComponent> defComponentMap;
    std::vector<defiNet> defNets;
  };

  static int defDieAreaCallback(defrCallbackType_e typ, defiBox *box, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defDieArea = *box;
    return 0;
  }

  static int defUnitCallback(defrCallbackType_e typ, double unit, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defUnit = unit;
    return 0;
  }

  static int defRowCallback(defrCallbackType_e typ, defiRow *row, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defRows.push_back(*row);
    return 0;
  }

  static int defComponentCallback(defrCallbackType_e typ, defiComponent *component, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defComponentMap.emplace(component->id(), *component);
    return 0;
  }

  static int defNetCallback(defrCallbackType_e typ, defiNet *net, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    // ignore io pin
    for (int i = 0; i < net->numConnections(); i++)
    {
      if (strcmp(net->instance(i), "PIN") == 0)
        return 0;
    }
    db->defNets.push_back(*net);
    return 0;
  }

  void parseDef(const std::string &defFilename, DEF_DATABASE *defdb)
  {
    // init def reader and set callback
    defrInit();
    defrSetDieAreaCbk(defDieAreaCallback);
    defrSetUnitsCbk(defUnitCallback);
    defrSetRowCbk(defRowCallback);
    defrSetComponentCbk(defComponentCallback);
    defrSetNetCbk(defNetCallback);

    // open the def file
    FILE *defFile = nullptr;
    if ((defFile = fopen(defFilename.c_str(), "r")) == 0)
    {
      LOG_CRITICAL("Couldn't open DEF File `{}`", defFilename);
      exit(1);
    }

    // invoke the def reader
    int res = defrRead(defFile, defFilename.c_str(), defdb, 1);
    fclose(defFile);
    if (res != 0)
    {
      LOG_CRITICAL("An error occurs when parsing DEF file `{}`", defFilename);
      exit(1);
    }
  }

  std::pair<int, int> getOrientSize(int orient, int w, int h)
  {
    switch (orient)
    {
    case 0:
    case 2:
    case 4:
    case 6:
      return std::make_pair(w, h);
    case 1:
    case 3:
    case 5:
    case 7:
      return std::make_pair(h, w);
    };

    LOG_WARN("Weird Orientation. Set North...");
    return std::make_pair(w, h);
  }

  std::pair<int, int> getOrientPoint(int orient, int x, int y, int w, int h)
  {
    switch (orient)
    {
    case 0: // North
      return std::make_pair(x, y);
    case 1: // West
      return std::make_pair(h-y, x);
    case 2: // South
      return std::make_pair(w-x, h-y);
    case 3: // Ease
      return std::make_pair(y, w-x);
    case 4: // Flipped North
      return std::make_pair(w-x, y);
    case 5: // Flipped West
      return std::make_pair(y, x);
    case 6: // Flipped South
      return std::make_pair(x, h-y);
    case 7: // Flipped East
      return std::make_pair(h-y, w-x);
    };

    LOG_WARN("Weird Orientation. Set North...");
    return std::make_pair(x, y);
  }

  std::shared_ptr<PlacerBase> Parser::lefdefToPlacerBase(const std::string &lefFilename, const std::string &defFilename)
  {
    auto pb = std::make_shared<PlacerBase>();

    auto tech = lefToTechnology(lefFilename);
    DEF_DATABASE defdb;
    LOG_TRACE("Parse DEF Begin");
    parseDef(defFilename, &defdb);
    LOG_TRACE("Parse DEF End");

    pb->techStor_.push_back(tech);
    pb->techs_.push_back(tech.get());
    pb->techNameMap_.emplace(tech->name(), tech.get());
    pb->instStor_.reserve(defdb.defComponentMap.size());
    pb->netStor_.reserve(defdb.defNets.size());
    size_t numPins = 0;
    for(const auto& net: defdb.defNets)
    {
      numPins += net.numConnections();
    }
    pb->pinStor_.reserve(numPins);

    LOG_TRACE("Process def site");

    LOG_TRACE("Process die and rows");
    // Process die and rows
    pb->dieStor_.emplace_back();
    pb->dieStor_.back().setName("def");
    pb->dieStor_.back().setDieBox(defdb.defDieArea.xl(), defdb.defDieArea.yl(), defdb.defDieArea.xh(), defdb.defDieArea.yh());
    pb->dieStor_.back().setMaxUtil(1.0f);
    pb->dieStor_.back().setTech(tech.get());
    // assume all row has the same height and width

    int rowStartX = INT_MAX, rowStartY = INT_MAX;
    for (const auto &defRow : defdb.defRows)
    {
      int lx = static_cast<int>(defRow.x());
      int ly = static_cast<int>(defRow.y());
      rowStartX = std::min(lx, rowStartX);
      rowStartY = std::min(ly, rowStartY);
    }
    int rowWidth = static_cast<int>(defdb.defRows.front().xNum() * tech->siteSizeX());
    int rowHeight = tech->siteSizeY();
    int rowRepeatCount = static_cast<int>(defdb.defRows.size());
    pb->dieStor_.back().setRowParams(rowStartX, rowStartY, rowWidth, rowHeight, rowRepeatCount);
    pb->dies_.push_back(&pb->dieStor_.back());
    pb->dieNameMap_.emplace("def", &pb->dieStor_.back());

    LOG_TRACE("Process def component");
    // Process def component
    for (const auto &compIt : defdb.defComponentMap)
    {
      int instLx = compIt.second.placementX();
      int instLy = compIt.second.placementY();
      const LibCell* libcell = tech->libCell(compIt.second.name());
      int orient = compIt.second.placementOrient();
      int w = libcell->sizeX();
      int h = libcell->sizeY();
      std::pair<int, int> instSize = getOrientSize(orient, w, h);

      pb->instStor_.emplace_back();
      pb->instStor_.back().setBox(instLx, instLy, instLx + instSize.first, instLy + instSize.second);
      pb->instStor_.back().setFixed(compIt.second.isFixed());
      pb->instStor_.back().setMacro(libcell->isMacro());
      pb->instStor_.back().setName(compIt.first);
      pb->instStor_.back().setLibCellName(libcell->name());
      pb->instNameMap_.emplace(compIt.first, &pb->instStor_.back());
    }

    LOG_TRACE("Process def net");
    // Process def net
    for(const auto& defNet : defdb.defNets)
    {
      pb->netStor_.emplace_back();
      for(int i = 0; i < defNet.numConnections(); i++)
      {
        auto instName = std::string(defNet.instance(i));
        auto pinName = std::string(defNet.pin(i));

        // find component
        const auto& compIt = defdb.defComponentMap.find(instName);

        // find libcell
        const LibCell* libcell = tech->libCell(compIt->second.name());
        int orient = compIt->second.placementOrient();

        // find pinOffset
        const LibPin* libpin = libcell->libPin(pinName);
        auto pinLoc = getOrientPoint(orient, libpin->x(), libpin->y(), libcell->sizeX(), libcell->sizeY());

        // find instance
        Instance* inst = pb->inst(instName);

        // add pin
        pb->pinStor_.emplace_back();
        pb->pinStor_.back().setName(pinName);
        pb->pinStor_.back().setInstance(inst);
        pb->pinStor_.back().setNet(&pb->netStor_.back());
        pb->pinStor_.back().updateLocation(inst, pinLoc.first, pinLoc.second);

        // assign inst & net
        pb->netStor_.back().addPin(&pb->pinStor_.back());
        inst->addPin(&pb->pinStor_.back());
      }
      pb->netStor_.back().updateBox();
    }

    // pb->placeMacrosArea_ = 0;
    // pb->placeStdcellsArea_ = 0;
    // pb->fixedMacrosArea_ = 0;
    // pb->fixedStdcellsArea_ = 0;

    // record pointers of pin, net, insts to vector<T *>
    for(auto& inst : pb->instStor_)
    {
      int64_t instArea = static_cast<int64_t>(inst.dx()) * static_cast<int64_t>(inst.dy());

      pb->insts_.push_back(&inst);
      pb->dieStor_.back().addInstance(&inst);
      if(inst.isFixed())
      {
        pb->fixedInsts_.push_back(&inst);
      }
      else
      {
        pb->placeInsts_.push_back(&inst);
      }
    }

    for(auto& pin : pb->pinStor_)
    {
      pb->pins_.push_back(&pin);
    }

    for(auto& net : pb->netStor_)
    {
      pb->nets_.push_back(&net);
    }

    return pb;
  }
}