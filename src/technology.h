#ifndef __REPLACE_TECHNOLOGY__
#define __REPLACE_TECHNOLOGY__

#include <string>
#include <unordered_map>

namespace replace
{
  class LibPin
  {
  public:
    LibPin() : x_(0), y_(0) {}
    LibPin(int x, int y) : x_(x), y_(y) {}
    ~LibPin() = default;

    void setLocation(int x, int y) { x_ = x; y_ = y; }
    int x() const { return x_; }
    int y() const { return y_; }

  private:
    int x_, y_;
  };

  class LibCell
  {
  public:
    LibCell() : sizeX_(0), sizeY_(0), isMacro_(false) {}
    LibCell(int w, int h, bool macro) : sizeX_(w), sizeY_(h), isMacro_(macro) {}
    ~LibCell() = default;

    void setMacro(bool on) { isMacro_ = on; }
    void setSize(int w, int h) { sizeX_ = w; sizeY_ = h; }
    bool isMacro() const { return isMacro_; }
    int sizeX() const { return sizeX_; }
    int sizeY() const { return sizeY_; }

    void addPin(const std::string& pinName, LibPin* pin) { pinNameMap_.emplace(pinName, pin); }
    LibPin* pin(const std::string& pinName) { return pinNameMap_.at(pinName); }
    const LibPin* pin(const std::string& pinName) const { return pinNameMap_.at(pinName); }
    size_t numPins() const { return pinNameMap_.size(); }
    int64_t area() const { return (int64_t)sizeX_ * sizeY_; }

  private:
    int sizeX_, sizeY_;
    bool isMacro_;
    std::unordered_map<std::string, LibPin*> pinNameMap_;
  };

  class Technology
  {
    friend class Parser;
  public:
    Technology() = default;
    ~Technology() = default;

    const std::vector<LibCell*>& libCells() const { return cells_; }
    const std::vector<LibCell*>& libStdCells() const { return stdCells_; }
    const std::vector<LibCell*>& libMacros() const { return macros_; }
    int siteSizeX() const { return siteSizeX_; }
    int siteSizeY() const { return siteSizeY_; }

    const LibCell* cell(const std::string& name) { return cellNameMap_.at(name); }

    void printInfo() const;

  private:
    std::vector<std::vector<LibPin>> pinStor_;
    std::vector<LibCell> cellStor_;

    std::vector<LibCell*> cells_;
    std::vector<LibCell*> stdCells_;
    std::vector<LibCell*> macros_;
    std::unordered_map<std::string, LibCell*> cellNameMap_;

    int siteSizeX_, siteSizeY_;
  };
}

#endif