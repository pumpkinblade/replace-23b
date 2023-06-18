#ifndef __REPLACE_PARSER__
#define __REPLACE_PARSER__

#include <string>
#include <memory>
#include "placerBase.h"

namespace replace
{
  class Parser
  {
  public:
    static std::shared_ptr<PlacerBase> FromLefDef(const std::string& lefFilename, const std::string& defFilename);
  };
}
#endif