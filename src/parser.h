#ifndef __REPLACE_PARSER__
#define __REPLACE_PARSER__

#include <string>
#include <memory>
#include "technology.h"
#include "placerBase.h"

namespace replace
{
  class Parser
  {
  public:
    static std::shared_ptr<Technology> LefToTechnology(const std::string& lefFilename);
    static std::shared_ptr<PlacerBase> LefDefToPlacerBase(const std::string& lefFilename, const std::string& defFilename);
  };
}
#endif