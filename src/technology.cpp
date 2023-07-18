#include "technology.h"
#include "log.h"

namespace replace
{
  /////////////////////////////////////////////////////////
  // LibPin

  LibPin::LibPin()
      : name_(), x_(0), y_(0)
  {
  }

  LibPin::LibPin(const std::string& name, int x, int y)
      : name_(name), x_(x), y_(y)
  {
  }

  ///////////////////////////////////////////////////////
  // LibCell

  LibCell::LibCell()
      : sizeX_(0), sizeY_(0), isMacro_(false)
  {
  }

  LibCell::LibCell(const std::string& name, int sizeX, int sizeY, bool isMacro)
      : name_(name), sizeX_(sizeX), sizeY_(sizeY), isMacro_(isMacro)
  {
  }

  LibPin* LibCell::libPin(const std::string& name) const
  {
    auto it = pinNameMap_.find(name);
    return it == pinNameMap_.end() ? nullptr : it->second;
  }

  void LibCell::addLibPin(LibPin* pin)
  {
    pins_.push_back(pin);
    pinNameMap_.emplace(pin->name(), pin);
  }

  /////////////////////////////////////////////////////////
  // Technology

  Technology::Technology()
      : siteSizeX_(0), siteSizeY_(0)
  {
  }

  Technology::Technology(const std::string& name)
      : name_(name), siteSizeX_(0), siteSizeY_(0)
  {
  }

  LibCell* Technology::libCell(const std::string& name) const
  {
    auto it = cellNameMap_.find(name);
    return it == cellNameMap_.end() ? nullptr : it->second;
  }

  void Technology::printInfo() const
  {
    int numLibStdCell = 0, numLibMacro = 0;
    int64_t maxStdArea = INT64_MIN;
    int64_t maxMacroArea = INT64_MIN;
    int64_t sumStdArea = 0;
    int64_t sumMacroArea = 0;

    for(const LibCell* cell : cells_)
    {
      int64_t area = (int64_t)cell->sizeX() * cell->sizeY();
      if(cell->isMacro())
      {
        maxMacroArea = std::max(area, maxMacroArea);
        sumMacroArea += area;
        numLibMacro++;
      }
      else
      {
        maxStdArea = std::max(area, maxStdArea);
        sumStdArea += area;
        numLibStdCell++;
      }
    }

    LOG_INFO("NumLibCells: {}", cells_.size());
    LOG_INFO("NumLibStdCells: {}", numLibStdCell);
    LOG_INFO("NumLibMacros: {}", numLibMacro);
    LOG_INFO("MaxStdArea: {}", maxStdArea);
    LOG_INFO("AvgStdArea: {}", sumStdArea / (double)numLibStdCell);
    LOG_INFO("MaxMacroArea: {}", maxMacroArea);
    LOG_INFO("AvgMacroArea: {}", sumMacroArea / (double)numLibMacro);

    int maxStdPin = INT_MIN;
    int maxMacroPin = INT_MIN;
    int sumStdPin = 0;
    int sumMacroPin = 0;
    for(const LibCell* cell : cells_)
    {
      if(cell->isMacro())
      {
        maxMacroPin = std::max((int)cell->libPins().size(), maxStdPin);
        sumMacroPin += (int)cell->libPins().size();
      }
      else
      {
        maxStdPin = std::max((int)cell->libPins().size(), maxMacroPin);
        sumStdPin += (int)cell->libPins().size();
      }
    }
    LOG_INFO("MaxStdPin: {}", maxStdPin);
    LOG_INFO("AvgStdPin: {}", sumStdPin / numLibStdCell);
    LOG_INFO("MaxMacroPin: {}", maxMacroPin);
    LOG_INFO("AvgMacroPin: {}", sumMacroPin / numLibMacro);
  }

}