#ifndef __REPLACE_PARSER__
#define __REPLACE_PARSER__

#include <string>
#include <memory>

namespace replace
{
  class Technology;
  class PlacerBase;

  class Parser
  {
  public:
    static std::shared_ptr<PlacerBase> txtToPlacerBase(const std::string& txtFilename);
  };
}
#endif
