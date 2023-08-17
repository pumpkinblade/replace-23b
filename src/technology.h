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
    friend class Parser;
  public:
    LibPin();
    LibPin(int id, int x, int y);
    ~LibPin() = default;

    int id() const { return id_; }
    void setId(int id) { id_ = id; }

    void setLocation(int x, int y) { x_ = x; y_ = y; }
    int x() const { return x_; }
    int y() const { return y_; }

  private:
    int id_;
    int x_, y_;
  };

  // LibCell can be macro or not, has shape, and named pins
  class LibCell
  {
    friend class Parser;
  public:
    LibCell();
    LibCell(int id, int sizeX, int sizeY, bool isMacro = false);
    ~LibCell() = default;

    int id() const { return id_; }
    void setId(int id) { id_ = id; }

    void setMacro(bool on) { isMacro_ = on; }
    bool isMacro() const { return isMacro_; }

    void setSize(int w, int h) { sizeX_ = w; sizeY_ = h; }
    int sizeX() const { return sizeX_; }
    int sizeY() const { return sizeY_; }

    const std::vector<LibPin*>& libPins() const { return pins_; }
    LibPin* libPin(int id) const { return pins_.at(id); }

  private:
    int id_;
    int sizeX_, sizeY_;
    bool isMacro_;
    std::vector<LibPin*> pins_;
  };

  // Technology is primarily a cellname-to-cellclass map
  class Technology
  {
    friend class Parser;
  public:
    Technology() = default;
    Technology(const std::string& name);
    ~Technology() = default;

    const std::string& name() const { return name_; }
    void setName(const std::string& name) { name_ = name; }

    const std::vector<LibCell*>& libCells() const { return cells_; }
    LibCell* libCell(int id) const { return cells_.at(id); }

    void printDebugInfo() const;

  private:
    std::string name_;
    std::vector<LibPin> pinStor_;
    std::vector<LibCell> cellStor_;
    std::vector<LibCell*> cells_;
  };
}

#endif