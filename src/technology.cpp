#include "technology.h"
#include "log.h"

namespace replace
{
  /////////////////////////////////////////////////////////
  // LibPin

  LibPin::LibPin()
      : id_(0), x_(0), y_(0)
  {
  }

  LibPin::LibPin(int id, int x, int y)
      : id_(id), x_(x), y_(y)
  {
  }

  ///////////////////////////////////////////////////////
  // LibCell

  LibCell::LibCell()
      : sizeX_(0), sizeY_(0), isMacro_(false)
  {
  }

  LibCell::LibCell(int id, int sizeX, int sizeY, bool isMacro)
      : id_(id), sizeX_(sizeX), sizeY_(sizeY), isMacro_(isMacro)
  {
  }
  
  /////////////////////////////////////////////////////////
  // Technology

  Technology::Technology(const std::string& name)
      : name_(name)
  {
  }

  void Technology::printDebugInfo() const
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

    LOG_DEBUG("NumLibCells: {}", cells_.size());
    LOG_DEBUG("NumLibStdCells: {}", numLibStdCell);
    LOG_DEBUG("NumLibMacros: {}", numLibMacro);
    LOG_DEBUG("MaxStdArea: {}", maxStdArea);
    LOG_DEBUG("AvgStdArea: {}", sumStdArea / (double)numLibStdCell);
    LOG_DEBUG("MaxMacroArea: {}", maxMacroArea);
    LOG_DEBUG("AvgMacroArea: {}", sumMacroArea / (double)numLibMacro);

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
    LOG_DEBUG("MaxStdPin: {}", maxStdPin);
    LOG_DEBUG("AvgStdPin: {}", sumStdPin / numLibStdCell);
    LOG_DEBUG("MaxMacroPin: {}", maxMacroPin);
    if(numLibMacro != 0){
      LOG_DEBUG("AvgMacroPin: {}", sumMacroPin / numLibMacro);
    }
  }

}