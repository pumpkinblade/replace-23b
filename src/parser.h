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
    static std::shared_ptr<PlacerBase> lefdefToPlacerBase(const std::string& lefFilename, const std::string& defFilename);
    // static std::shared_ptr<Placer23b> txtToPlacer23b(const std::string& txtFilename);
    static std::shared_ptr<PlacerBase> txtToPlacerBase(const std::string& txtFilename);
  private:
    static std::shared_ptr<Technology> lefToTechnology(const std::string& lefFilename);
  };
}
#endif