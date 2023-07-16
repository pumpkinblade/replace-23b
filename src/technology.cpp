#include "technology.h"
#include "log.h"

namespace replace
{
  /////////////////////////////////////////////////////////
  // Technology

  void Technology::printInfo() const
  {
    LOG_INFO("NumLibCells: {}", cells_.size());
    LOG_INFO("NumLibStdCells: {}", stdCells_.size());
    LOG_INFO("NumLibMacros: {}", macros_.size());

    int64_t maxStdArea = INT64_MIN;
    int64_t maxMacroArea = INT64_MIN;
    int64_t sumStdArea = 0;
    int64_t sumMacroArea = 0;
    for(const LibCell* cell : cells_)
    {
      int64_t area = cell->area();
      if(cell->isMacro())
      {
        maxMacroArea = std::max(area, maxMacroArea);
        sumMacroArea += area;
      }
      else
      {
        maxStdArea = std::max(area, maxStdArea);
        sumStdArea += area;
      }
    }
    LOG_INFO("MaxStdArea: {}", maxStdArea);
    LOG_INFO("AvgStdArea: {}", sumStdArea / (double)stdCells_.size());
    LOG_INFO("MaxMacroArea: {}", maxMacroArea);
    LOG_INFO("AvgMacroArea: {}", sumMacroArea / (double)macros_.size());

    int maxStdPin = INT_MIN;
    int maxMacroPin = INT_MIN;
    int sumStdPin = 0;
    int sumMacroPin = 0;
    for(const LibCell* cell : cells_)
    {
      if(cell->isMacro())
      {
        maxMacroPin = std::max((int)cell->numPins(), maxStdPin);
        sumMacroPin += (int)cell->numPins();
      }
      else
      {
        maxStdPin = std::max((int)cell->numPins(), maxMacroPin);
        sumStdPin += (int)cell->numPins();
      }
    }
    LOG_INFO("MaxStdPin: {}", maxStdPin);
    LOG_INFO("AvgStdPin: {}", sumStdPin / (float)stdCells_.size());
    LOG_INFO("MaxMacroPin: {}", maxMacroPin);
    LOG_INFO("AvgMacroPin: {}", sumMacroPin / (float)macros_.size());
  }

}