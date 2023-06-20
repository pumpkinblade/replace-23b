#include "parser.h"
#include "log.h"
#include <algorithm>
#include <lefrReader.hpp>
#include <defrReader.hpp>

namespace replace
{
  struct LEF_DATABASE
  {
    lefiUnits lefUnit;
    std::unordered_map<std::string, lefiSite> lefSiteMap;
    std::unordered_map<std::string, lefiMacro> lefMacroMap;
    std::unordered_map<std::string, int> lefMacroPinIdxMap;
    std::vector<std::vector<lefiPin>> lefMacroPins;
  };

  static int LefUnitCallback(lefrCallbackType_e typ, lefiUnits *unit, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefUnit = *unit;
    return 0;
  }

  static int LefSiteCallback(lefrCallbackType_e typ, lefiSite *site, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefSiteMap.emplace(site->name(), *site);
    return 0;
  }

  static int LefMacroCallback(lefrCallbackType_e typ, lefiMacro *macro, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefMacroMap.emplace(macro->name(), *macro);
    return 0;
  }

  static int LefMacroBeginCallback(lefrCallbackType_e typ, const char *name, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefMacroPinIdxMap.emplace(name, db->lefMacroPins.size());
    db->lefMacroPins.emplace_back();
    return 0;
  }

  static int LefPinCallback(lefrCallbackType_e typ, lefiPin *pin, lefiUserData data)
  {
    LEF_DATABASE *db = reinterpret_cast<LEF_DATABASE *>(data);
    db->lefMacroPins.back().push_back(*pin);
    return 0;
  }

  void ParseLef(const std::string &lefFilename, LEF_DATABASE *lefdb)
  {
    // init lef reader and set callback
    lefrInit();
    lefrSetUnitsCbk(LefUnitCallback);
    lefrSetSiteCbk(LefSiteCallback);
    lefrSetMacroBeginCbk(LefMacroBeginCallback);
    lefrSetMacroCbk(LefMacroCallback);
    lefrSetPinCbk(LefPinCallback);

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

  struct DEF_DATABASE
  {
    defiBox defDieArea;
    double defUnit;
    std::vector<defiRow> defRows;
    std::unordered_map<std::string, defiComponent> defComponentMap;
    std::vector<defiNet> defNets;
  };

  static int DefDieAreaCallback(defrCallbackType_e typ, defiBox *box, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defDieArea = *box;
    return 0;
  }

  static int DefUnitCallback(defrCallbackType_e typ, double unit, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defUnit = unit;
    return 0;
  }

  static int DefRowCallback(defrCallbackType_e typ, defiRow *row, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defRows.push_back(*row);
    return 0;
  }

  static int DefComponentCallback(defrCallbackType_e typ, defiComponent *component, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defComponentMap.emplace(component->id(), *component);
    return 0;
  }

  static int DefNetCallback(defrCallbackType_e typ, defiNet *net, defiUserData data)
  {
    DEF_DATABASE *db = reinterpret_cast<DEF_DATABASE *>(data);
    db->defNets.push_back(*net);
    return 0;
  }

  void ParseDef(const std::string &defFilename, DEF_DATABASE *defdb)
  {
    // init def reader and set callback
    defrInit();
    defrSetDieAreaCbk(DefDieAreaCallback);
    defrSetUnitsCbk(DefUnitCallback);
    defrSetRowCbk(DefRowCallback);
    defrSetComponentCbk(DefComponentCallback);
    defrSetNetCbk(DefNetCallback);

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

  std::pair<double, double> GetPinOffset(const lefiPin &pin)
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

  std::pair<int, int> GetOrientSize(int orient, int w, int h)
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

  std::pair<int, int> GetOrientPoint(int orient, int x, int y, int w, int h)
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

  std::shared_ptr<PlacerBase> Parser::FromLefDef(const std::string &lefFilename, const std::string &defFilename)
  {
    LEF_DATABASE lefdb;
    DEF_DATABASE defdb;

    LOG_TRACE("Parse LEF Begin");
    ParseLef(lefFilename, &lefdb);
    LOG_TRACE("Parse LEF End");
    LOG_TRACE("Parse DEF Begin");
    ParseDef(defFilename, &defdb);
    LOG_TRACE("Parse DEF End");

    LOG_INFO("Unit: {} {}", lefdb.lefUnit.databaseName(), lefdb.lefUnit.databaseNumber());
    LOG_INFO("Num Site: {}", lefdb.lefSiteMap.size());
    LOG_INFO("Num Macro: {}", lefdb.lefMacroMap.size());

    LOG_INFO("Die Area: [{}, {}] x [{}, {}]", defdb.defDieArea.xl(), defdb.defDieArea.xh(), defdb.defDieArea.yl(), defdb.defDieArea.yh());
    LOG_INFO("DEF Unit: {}", defdb.defUnit);
    LOG_INFO("Num Row: {}", defdb.defRows.size());
    LOG_INFO("Num Component: {}", defdb.defComponentMap.size());
    LOG_INFO("Num Net: {}", defdb.defNets.size());

    auto pb = std::make_shared<PlacerBase>();
    pb->instStor_.reserve(defdb.defComponentMap.size());
    pb->netStor_.reserve(defdb.defNets.size());
    int numPins = 0;
    for(const auto& net: defdb.defNets)
    {
      numPins += net.numConnections();
    }
    pb->pinStor_.reserve(numPins);

    LOG_TRACE("Process Site");
    // assume that every row uses same site
    const char* siteName = defdb.defRows.front().macro();
    const auto& site = lefdb.lefSiteMap[siteName];
    int siteWidth = static_cast<int>(site.sizeX() * defdb.defUnit);
    int siteHeight = static_cast<int>(site.sizeY() * defdb.defUnit);
    pb->siteSizeX_ = siteWidth;
    pb->siteSizeY_ = siteHeight;

    LOG_TRACE("Process die and rows");
    // Process die and rows
    pb->die_.setDieBox(defdb.defDieArea.xl(), defdb.defDieArea.yl(), defdb.defDieArea.xh(), defdb.defDieArea.yh());
    for (const auto &defRow : defdb.defRows)
    {
      int lx = static_cast<int>(defRow.x());
      int ly = static_cast<int>(defRow.y());
      int ux = lx + static_cast<int>(defRow.xNum()) * siteWidth;
      int uy = ly + siteHeight;
      pb->die_.addRow(Row(siteWidth, lx, ly, ux, uy));
    }
    pb->die_.updateCoreBox();

    LOG_TRACE("Process lef macro");
    // Process lef macro
    std::unordered_map<std::string, std::pair<int, int>> pinOffsets; // $macroName-$pinName -> (pinOffsetX, pinOffsetY)
    for (const auto &lefMacroIt : lefdb.lefMacroMap)
    {
      int idx = lefdb.lefMacroPinIdxMap[lefMacroIt.first];
      for (const auto &pin : lefdb.lefMacroPins[idx])
      {
        std::string key = lefMacroIt.first + std::string("-") + pin.name();
        std::pair<double, double> offset = GetPinOffset(pin);
        int offsetX = static_cast<int>(offset.first * defdb.defUnit);
        int offsetY = static_cast<int>(offset.second * defdb.defUnit);
        pinOffsets.emplace(key, std::make_pair(offsetX, offsetY));
      }
    }

    LOG_TRACE("Process def component");
    // Process def component
    std::unordered_map<std::string, int> instExtIds;
    for (const auto &compIt : defdb.defComponentMap)
    {
      int instLx = compIt.second.placementX();
      int instLy = compIt.second.placementY();
      const auto& macro = lefdb.lefMacroMap[compIt.second.name()];
      int orient = compIt.second.placementOrient();
      int w = static_cast<int>(macro.sizeX() * defdb.defUnit);
      int h = static_cast<int>(macro.sizeY() * defdb.defUnit);
      std::pair<int, int> instSize = GetOrientSize(orient, w, h);

      Instance inst;
      inst.setBox(instLx, instLy, instLx + instSize.first, instLy + instSize.second);
      inst.setFixed(compIt.second.isFixed());
      inst.setDummy(false);
      inst.setExtId(static_cast<int>(pb->instStor_.size()));
      pb->instStor_.push_back(inst);

      instExtIds.emplace(compIt.second.id(), inst.extId());
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

        // find macro
        auto macroName = std::string(compIt->second.name());
        const auto& macro = lefdb.lefMacroMap[macroName];
        int orient = compIt->second.placementOrient();

        // find pinOffset
        std::pair<int, int> pinOffset = pinOffsets[macroName + "-" + pinName];
        int x = pinOffset.first;
        int y = pinOffset.second;
        int w = static_cast<int>(macro.sizeX() * defdb.defUnit);
        int h = static_cast<int>(macro.sizeY() * defdb.defUnit);
        pinOffset = GetOrientPoint(orient, x, y, w, h);

        // find instance extId
        int extId = instExtIds[instName];

        // add pin
        pb->pinStor_.emplace_back();
        pb->pinStor_.back().setInstance(&pb->instStor_[extId]);
        pb->pinStor_.back().setNet(&pb->netStor_.back());
        pb->pinStor_.back().updateLocation(&pb->instStor_[extId], pinOffset.first, pinOffset.second);
        pb->netStor_.back().addPin(&pb->pinStor_.back());
        pb->instStor_[extId].addPin(&pb->pinStor_.back());
      }
      pb->netStor_.back().updateBox();
    }

    pb->placeInstsArea_ = 0;
    pb->macroInstsArea_ = 0;
    pb->stdInstsArea_ = 0;
    pb->nonPlaceInstsArea_ = 0;
    for(auto& inst : pb->instStor_)
    {
      int64_t instArea = static_cast<int64_t>(inst.dx()) * static_cast<int64_t>(inst.dy());

      pb->insts_.push_back(&inst);
      if(inst.isFixed())
      {
        pb->fixedInsts_.push_back(&inst);
        pb->nonPlaceInsts_.push_back(&inst);
        pb->nonPlaceInstsArea_ += instArea;
      }
      else if(inst.isDummy())
      {
        pb->dummyInsts_.push_back(&inst);
        pb->nonPlaceInsts_.push_back(&inst);
        pb->nonPlaceInstsArea_ += instArea;
      }
      else
      {
        pb->placeInsts_.push_back(&inst);
        pb->placeInstsArea_ += instArea;
      }
      
      if(inst.dy() > pb->siteSizeY_ * 6)
        pb->macroInstsArea_ += instArea;
      else
        pb->stdInstsArea_ += instArea;
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