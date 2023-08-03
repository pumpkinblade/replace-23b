#ifndef __REPLACE_TECHNOLOGY__
#define __REPLACE_TECHNOLOGY__

#include <string>
#include <vector>
#include <unordered_map>
#include <vector>

namespace replace
{
  class LibPin
  {
  public:
    LibPin();
    LibPin(const std::string& name, int x, int y);
    ~LibPin() = default;

    const std::string& name() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    void setLocation(int x, int y) { x_ = x; y_ = y; }
    int x() const { return x_; }
    int y() const { return y_; }

  private:
    std::string name_;
    int x_, y_;
  };

  // LibCell can be macro or not, has shape, and named pins
  class LibCell
  {
  public:
    LibCell();
    LibCell(const std::string& name, int sizeX, int sizeY, bool isMacro = false);
    ~LibCell() = default;

    const std::string& name() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    void setMacro(bool on) { isMacro_ = on; }
    bool isMacro() const { return isMacro_; }

    void setSize(int w, int h) { sizeX_ = w; sizeY_ = h; }
    int sizeX() const { return sizeX_; }
    int sizeY() const { return sizeY_; }

    const std::vector<LibPin*>& libPins() const { return pins_; }
    LibPin* libPin(const std::string& name) const;
    void addLibPin(LibPin* pin);

  private:
    std::string name_;
    int sizeX_, sizeY_;
    bool isMacro_;
    std::vector<LibPin*> pins_;
    std::unordered_map<std::string, LibPin*> pinNameMap_;
  };

  // Technology is primarily a cellname-to-cellclass map
  class Technology
  {
    friend class Parser;
  public:
    Technology();
    Technology(const std::string& name);
    ~Technology() = default;

    const std::string& name() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    const std::vector<LibCell*>& libCells() const { return cells_; }
    LibCell* libCell(const std::string& name) const;

    int siteSizeX() const { return siteSizeX_; }
    int siteSizeY() const { return siteSizeY_; }

    void printDebugInfo() const;

  private:
    std::string name_;

    std::vector<LibPin> pinStor_;
    std::vector<LibCell> cellStor_;

    std::vector<LibCell*> cells_;
    std::unordered_map<std::string, LibCell*> cellNameMap_;

    // For lef/def parser
    int siteSizeX_, siteSizeY_;
  };
}

#endif