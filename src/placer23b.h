#ifndef __REPLACE_PLACER_23B__
#define __REPLACE_PLACER_23B__

#include "technology.h"
#include "placerBase.h"

namespace replace
{
  class Instance23b
  {
  public:
    Instance23b()
        : libCellName_() {}
    Instance23b(const std::string& cellName)
        : libCellName_(cellName) {}
    ~Instance23b() = default;

    const std::string& libCellName() const { return libCellName_; }

  private:
    std::string libCellName_;
  };

  class Pin23b
  {
  public:
    Pin23b()
        : instName_(), pinName_() {}
    Pin23b(const std::string& instName, const std::string& pinName)
        : instName_(instName), pinName_(pinName) {}
    ~Pin23b() = default;

    const std::string& instName() const { return instName_; }
    const std::string& pinName() const { return pinName_; }

  private:
    std::string instName_;
    std::string pinName_;
  };

  class Net23b
  {
  public:
    Net23b() = default;
    ~Net23b() = default;

    void addPin(Pin23b* pin) { pins_.push_back(pin); }

    const std::vector<Pin23b*>& pins() const { return pins_; }

  private:
    std::vector<Pin23b*> pins_;
  };

  class Placer23b
  {
    friend class Parser;
  public:
    const Die &topDie() const { return topDie_; }
    const Die &boottomDie() const { return bottomDie_; }
    float topDieMaxUtil() const { return topUtil_; }
    float bottomDieMaxUtil() const { return bottomUtil_; }
    Technology* topDieTechnology() const { return topTech_; }
    Technology* bottomDieTechnology() const { return bottomTech_; }

    int terminalSizeX() const { return termSizeX_; }
    int terminalSizeY() const { return termSizeY_; }
    int terminalSpacing() const { return termSpace_; }
    int terminalCost() const { return termCost_; }

    const std::vector<Instance23b*>& insts() const { return insts_; }
    const std::vector<Pin23b*>& pins() const { return pins_; }
    const std::vector<Net23b*>& nets() const { return nets_; }
    Instance23b* inst(const std::string& name) const { return instNameMap_.at(name); }
    Net23b* net(const std::string& name) const { return netNameMap_.at(name); }
    Technology* technology(const std::string& name) { return techNameMap_.at(name); }

    void printInfo() const;

  private:
    std::vector<Technology> techStor_;
    std::vector<Instance23b> instStor_;
    std::vector<Pin23b> pinStor_;
    std::vector<Net23b> netStor_;

    Die topDie_;
    Die bottomDie_;
    Technology* topTech_;
    Technology* bottomTech_;
    float topUtil_;
    float bottomUtil_;

    int termSizeX_;
    int termSizeY_;
    int termSpace_;
    int termCost_;

    std::vector<Instance23b*> insts_;
    std::vector<Pin23b*> pins_;
    std::vector<Net23b*> nets_;
    std::unordered_map<std::string, Instance23b*> instNameMap_;
    std::unordered_map<std::string, Net23b*> netNameMap_;
    std::unordered_map<std::string, Technology*> techNameMap_;
  };
}

#endif