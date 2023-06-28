#ifndef __REPLACE_PARSER__
#define __REPLACE_PARSER__

#include <string>
#include <memory>
#include "technology.h"
#include "placerBase.h"
#include "placer23b.h"

namespace replace
{
  class Parser
  {
  public:
    static std::shared_ptr<PlacerBase> LefDefToPlacerBase(const std::string& lefFilename, const std::string& defFilename);
    static std::shared_ptr<Placer23b> TxtToPlacer23b(const std::string& txtFilename);
  private:
    static std::shared_ptr<Technology> LefToTechnology(const std::string& lefFilename);
  };
}
#endif