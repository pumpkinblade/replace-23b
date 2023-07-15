#ifndef __REPLACE_PLOT__
#define __REPLACE_PLOT__

#include <string>
#include <memory>

namespace replace
{
  class PlacerBase;
  class NesterovBase;

  class Plotter;

  class PlotVars
  {
  public:
    int minLength;
    int xMargin;
    int yMargin;
  };

  enum class PlotNesterovType
  {
    GCell,
    Bin,
    Arrow
  };

  class Plot
  {
  public:
    static void init(const PlotVars &vars);

    static void plot(const PlacerBase *pb,
                     const std::string &imgDir,
                     const std::string &prefix);
    static void plot(const NesterovBase *nb,
                     PlotNesterovType type,
                     const std::string &imgDir,
                     const std::string &prefix);

  private:
    static std::unique_ptr<Plotter> splotter_;
  };
}

#endif
