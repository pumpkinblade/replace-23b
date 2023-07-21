#ifndef __REPLACE_OUTPUT_WRITER__
#define __REPLACE_OUTPUT_WRITER__

#include <string>

namespace replace
{
  class PlacerBase;
  class Instance;

  class OutputWriter
  {
  public:
    static void write(PlacerBase* pb, std::string& filename);
  };
}

#endif