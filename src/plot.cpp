#include "plot.h"
#include "placerBase.h"
#include "nesterovBase.h"
#include "log.h"
#include <array>
#include <vector>
#include <CImg.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace replace
{
  using Color = std::array<unsigned char, 3>;
  constexpr Color g_red = {255, 0, 0};
  constexpr Color g_green = {0, 255, 0};
  constexpr Color g_blue = {0, 0, 255};
  constexpr Color g_pink = {255, 0, 255};
  constexpr Color g_cyron = {0, 255, 255};
  constexpr Color g_yellow = {255, 255, 0};
  constexpr Color g_black = {0, 0, 0};
  constexpr Color g_white = {255, 255, 255};

  using namespace cimg_library;
  using CImgObj = CImg<unsigned char>;

  class Plotter
  {
  public:
    Plotter(const PlotVars &vars);
    ~Plotter() = default;

    void plot(const PlacerBase *pb,
              const std::string &imgDir,
              const std::string &prefix);
    void plot(const NesterovBase *nb,
              PlotNesterovType type,
              const std::string &imgDir,
              const std::string &prefix);

  private:
    void initContext(const Die *die);
    void initContext(const BinGrid *bg);

    int getImageX(int origX) const;
    int getImageY(int origY) const;

    void drawInstances(const Die *die, float opacity);
    void drawGCells(const BinGrid *bg, float opacity);
    void drawBins(const BinGrid *bg, float opacity);
    void drawArrows(const BinGrid *bg, float opacity);

    void cimgDrawArrow(int x1, int y1, int x3, int y3, float thick, const Color &color, float opacity);
    void cimgWriteJpeg(const std::string& name, unsigned int quality);
    void cimgWritePng(const std::string& name);

  private:
    int minLength_;
    int xMargin_;
    int yMargin_;

    int imageWidth_;
    int imageHeight_;
    int origLx_;
    int origLy_;
    int origWidth_;
    int origHeight_;
    float unitX_;
    float unitY_;

    std::unique_ptr<CImgObj> img_;
  };

  ///////////////////////////////////////////
  // Plotter

  Plotter::Plotter(const PlotVars &vars)
      : minLength_(0), xMargin_(0), yMargin_(0),
        imageWidth_(0), imageHeight_(0),
        origLx_(0), origLy_(0), origWidth_(0), origHeight_(0),
        unitX_(0), unitY_(0),
        img_(nullptr)
  {
    minLength_ = vars.minLength;
    xMargin_ = vars.xMargin;
    yMargin_ = vars.yMargin;
  }

  void Plotter::plot(const PlacerBase *pb,
                     const std::string &imgDir,
                     const std::string &prefix)
  {
    for (const Die* die : pb->dies())
    {
      initContext(die);

      drawInstances(die, 0.7f);

      std::string title = prefix + "_" + die->name();
      img_->draw_text(0, 0, title.c_str(), g_black.data(), NULL, 1, 30);

      std::string saveName = imgDir + "/" + title + ".jpg";
      // img_->save_jpeg(saveName.c_str(), 70);
      cimgWriteJpeg(saveName, 70);

      LOG_TRACE("{} has been saved.", saveName);
    }
  }

  void Plotter::plot(const NesterovBase *nb,
                     PlotNesterovType type,
                     const std::string &imgDir,
                     const std::string &prefix)
  {
    for (BinGrid *bg : nb->binGrids())
    {
      initContext(bg);

      switch (type)
      {
      case PlotNesterovType::GCell:
        drawGCells(bg, 0.7f);
        break;
      case PlotNesterovType::Bin:
        drawBins(bg, 0.7f);
        break;
      case PlotNesterovType::Arrow:
        drawBins(bg, 0.7f);
        drawArrows(bg, 0.7f);
        break;
      default:
        break;
      }

      std::string title = prefix + "_" + bg->die()->name();
      img_->draw_text(0, 0, title.c_str(), g_black.data(), NULL, 1, 30);

      std::string saveName = imgDir + "/" + title + ".jpg";
      // img_->save_jpeg(saveName.c_str(), 70);
      cimgWriteJpeg(saveName, 70);

      LOG_TRACE("{} has been saved.", saveName);
    }
  }

  void Plotter::initContext(const Die *die)
  {
    origLx_ = die->dieLx();
    origLy_ = die->dieLy();
    origWidth_ = die->dieDx();
    origHeight_ = die->dieDy();

    int newWidth, newHeight;
    // imageWidth & height setting
    // Set minimum length of picture as minLength
    if (origWidth_ < origHeight_)
    {
      newHeight = static_cast<int>(origHeight_ / (origWidth_ / (float)minLength_));
      newWidth = minLength_;
    }
    else
    {
      newWidth = static_cast<int>(origWidth_ / (origHeight_ / (float)minLength_));
      newHeight = minLength_;
    }

    if (img_ != nullptr && newWidth == imageWidth_ && newHeight == imageHeight_)
    {
      img_->fill(255);
    }
    else
    {
      CImgObj *img = new CImgObj(newWidth + 2 * xMargin_,
                                 newHeight + 2 * yMargin_,
                                 1, 3, 255);
      img_.reset(img);
      imageWidth_ = newWidth;
      imageHeight_ = newHeight;
    }

    // scaling
    unitX_ = imageWidth_ / (float)origWidth_;
    unitY_ = imageHeight_ / (float)origHeight_;
  }

  void Plotter::initContext(const BinGrid *bg)
  {
    initContext(bg->die());
  }

  int Plotter::getImageX(int origX) const
  {
    return static_cast<int>((origX - origLx_) * unitX_ + xMargin_);
  }

  int Plotter::getImageY(int origY) const
  {
    return static_cast<int>((origHeight_ - (origY - origLy_)) * unitY_ + yMargin_);
  }

  void Plotter::drawInstances(const Die *die, float opacity)
  {
    for (const Instance *inst : die->insts())
    {
      int x1 = getImageX(inst->lx());
      int y1 = getImageY(inst->ly());
      int x3 = getImageX(inst->ux());
      int y3 = getImageY(inst->uy());

      if(inst->isMacro())
      {
        img_->draw_rectangle(x1, y1, x3, y3, g_cyron.data(), opacity);
      }
      else if (inst->isFixed())
      {
        img_->draw_rectangle(x1, y1, x3, y3, g_blue.data(), opacity);
      }
      else
      {
        img_->draw_rectangle(x1, y1, x3, y3, g_red.data(), opacity);
      }
    }
  }

  void Plotter::drawGCells(const BinGrid *bg, float opacity)
  {
    for (const GCell *cell : bg->gCells())
    {
      int cx = cell->cx();
      int cy = cell->cy();
      int x1 = getImageX(cx - cell->dx() / 2);
      int y1 = getImageY(cy - cell->dy() / 2);
      int x2 = getImageX(cx + cell->dx() / 2);
      int y2 = getImageY(cy + cell->dy() / 2);

      if (cell->isInstance())
      {
        img_->draw_rectangle(x1, y1, x2, y2, g_blue.data(), opacity);
      }
    }
  }

  void Plotter::drawBins(const BinGrid *bg, float opacity)
  {
    for (const Bin *bin : bg->bins())
    {
      int x1 = getImageX(bin->lx());
      int y1 = getImageY(bin->ly());
      int x2 = getImageX(bin->ux());
      int y2 = getImageY(bin->uy());

      int color = static_cast<int>(bin->density() * 50 + 20);
      color = 255 - std::max(std::min(color, 255), 20);
      unsigned char drawColor = static_cast<unsigned char>(color);
      unsigned char c[3] = {drawColor, drawColor, drawColor};
      img_->draw_rectangle(x1, y1, x2, y2, c, opacity);
    }
  }

  void Plotter::drawArrows(const BinGrid *bg, float opacity)
  {
    int binMaxX = bg->binCntX();
    int binMaxY = bg->binCntY();
    int arrowSpacing = (binMaxX / 16 <= 0) ? 1 : binMaxX / 16;

    // below is essential for extracting e?Max
    float exMax = std::numeric_limits<float>::epsilon();
    float eyMax = std::numeric_limits<float>::epsilon();
    for (int i = 0; i < binMaxX; i += arrowSpacing)
    {
      for (int j = 0; j < binMaxY; j += arrowSpacing)
      {
        const Bin *bin = bg->bins()[binMaxX * j + i];
        float newEx = fabs(bin->electroForceX());
        float newEy = fabs(bin->electroForceY());

        exMax = (exMax < newEx) ? newEx : exMax;
        eyMax = (eyMax < newEy) ? newEy : eyMax;
      }
    }

    for (int i = 0; i < binMaxX; i += arrowSpacing)
    {
      for (int j = 0; j < binMaxY; j += arrowSpacing)
      {
        Bin *bin = bg->bins()[binMaxX * j + i];
        int signX = (bin->electroForceX() > 0) ? 1 : -1;
        int signY = (bin->electroForceY() > 0) ? 1 : -1;

        float newVx = std::fabs(bin->electroForceX());
        float newVy = std::fabs(bin->electroForceY());

        int x1 = bin->cx();
        int y1 = bin->cy();

        float dx = signX * newVx / exMax;
        float dy = signY * newVy / eyMax;

        //        float theta = atan(dy / dx);
        float length = 5.0f * std::sqrt(std::pow((float)bg->binSizeX(), 2.0f) +
                                        std::pow((float)bg->binSizeY(), 2.0f));

        int x3 = static_cast<int>(x1 + dx * length);
        int y3 = static_cast<int>(y1 + dy * length);

        int drawX1 = getImageX(x1);
        int drawY1 = getImageY(y1);

        int drawX3 = getImageX(x3);
        int drawY3 = getImageY(y3);

        // img_->draw_arrow(drawX1, drawY1, drawX3, drawY3, g_black.data(), opacity);
        cimgDrawArrow(drawX1, drawY1, drawX3, drawY3, 20.f, g_red, opacity);
      }
    }
  }

  void Plotter::cimgDrawArrow(int x1, int y1, int x3, int y3, float thick, const Color &color, float opacity)
  {
    // ARROW HEAD DRAWING
    float arrowHeadSize = static_cast<float>(thick);
    float theta = std::atan((y3 - y1) / (float)(x3 - x1));
    float invTheta = std::atan(-(x3 - x1) / (float)(y3 - y1));

    // ARROW RECT DRAWING
    int ldX = static_cast<int>(x1 - thick / 4.f * std::cos(invTheta));
    int ldY = static_cast<int>(y1 - thick / 4.f * std::sin(invTheta));
    int rdX = static_cast<int>(x1 + thick / 4.f * std::cos(invTheta));
    int rdY = static_cast<int>(y1 + thick / 4.f * std::sin(invTheta));

    int luX = static_cast<int>(x3 - thick / 4.f * std::cos(invTheta));
    int luY = static_cast<int>(y3 - thick / 4.f * std::sin(invTheta));
    int ruX = static_cast<int>(x3 + thick / 4.f * std::cos(invTheta));
    int ruY = static_cast<int>(y3 + thick / 4.f * std::sin(invTheta));

    cimg_library::CImg<int> rectPoints(4, 2);
    rectPoints(0, 0) = ldX;
    rectPoints(0, 1) = ldY;
    rectPoints(1, 0) = rdX;
    rectPoints(1, 1) = rdY;
    rectPoints(2, 0) = ruX;
    rectPoints(2, 1) = ruY;
    rectPoints(3, 0) = luX;
    rectPoints(3, 1) = luY;

    img_->draw_polygon(rectPoints, color.data());

    cimg_library::CImg<int> headPoints(3, 2);
    int lPointX = static_cast<int>(x3 - thick / 2.f * std::cos(invTheta));
    int lPointY = static_cast<int>(y3 - thick / 2.f * std::sin(invTheta));
    int rPointX = static_cast<int>(x3 + thick / 2.f * std::cos(invTheta));
    int rPointY = static_cast<int>(y3 + thick / 2.f * std::sin(invTheta));

    int uPointX = (x3 - x1) >= 0 ? static_cast<int>(x3 + thick * std::cos(theta))
                                 : static_cast<int>(x3 - thick * std::cos(theta));
    int uPointY = (x3 - x1) >= 0 ? static_cast<int>(y3 + thick * std::sin(theta))
                                 : static_cast<int>(y3 - thick * std::sin(theta));

    //    int uPointX = x3 + 1.0*thick*cos(theta);
    //    int uPointY = y3 + 1.0*thick*sin(theta);

    headPoints(0, 0) = lPointX;
    headPoints(0, 1) = lPointY;
    headPoints(1, 0) = rPointX;
    headPoints(1, 1) = rPointY;
    headPoints(2, 0) = uPointX;
    headPoints(2, 1) = uPointY;

    img_->draw_polygon(headPoints, color.data());
  }

  void Plotter::cimgWriteJpeg(const std::string& name, unsigned int quality)
  {
    int w = img_->width();
    int h = img_->height();
    img_->permute_axes("cxyz");
    stbi_write_jpg(name.c_str(), w, h, 3, img_->data(), 70);
    img_->permute_axes("yzcx");
  }

  void Plotter::cimgWritePng(const std::string& name)
  {
    int w = img_->width();
    int h = img_->height();
    img_->permute_axes("cxyz");
    stbi_write_png(name.c_str(), w, h, 3, img_->data(), 3 * w);
    img_->permute_axes("yzcx");
  }

  ////////////////////////////////////////////
  // Plot

  std::unique_ptr<Plotter> Plot::splotter_;

  void Plot::init(const PlotVars &vars)
  {
    Plotter *plotter = new Plotter(vars);
    splotter_.reset(plotter);
  }

  void Plot::plot(const PlacerBase *pb,
                  const std::string &imgDir,
                  const std::string &prefix)
  {
    if (splotter_ == nullptr)
    {
      LOG_WARN("Couldn't find suitable plotter. This process will be skipped.");
    }
    else
    {
      splotter_->plot(pb, imgDir, prefix);
    }
  }

  void Plot::plot(const NesterovBase *nb,
                  PlotNesterovType type,
                  const std::string &imgDir,
                  const std::string &prefix)
  {
    if (splotter_ == nullptr)
    {
      LOG_WARN("Couldn't find suitable plotter. This process will be skipped.");
    }
    else
    {
      splotter_->plot(nb, type, imgDir, prefix);
    }
  }
}
