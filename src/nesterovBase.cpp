#include "nesterovBase.h"
#include "placerBase.h"
#include "fft.h"
#include "log.h"

#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>

namespace replace
{
  static double getOverlapArea(const Bin *bin, const Instance *inst);
  static double getOverlapDensityArea(const Bin *bin, const GCell *cell);
  static double fastExp(double exp);
  static double bellShape(double theta_i, double theta);
  static double bellShape1(double theta_i);
  static double bellShape2(double theta_i);
  static double bellShapeDerive(double theta_i, double theta);
  static double bellShapeDerive1(double theta_i);
  static double bellShapeDerive2(double theta_i);

  ////////////////////////////////////////////////
  // GCell

  GCell::GCell()
      : inst_(nullptr), isMacro_(false),
        cx_(0), cy_(0), dx_(0), dy_(0), theta_(0),
        densityScale_(0)
  {
  }

  GCell::GCell(Instance *inst)
      : GCell()
  {
    setInstance(inst);
  }

  GCell::GCell(double cx, double cy, double dx, double dy)
      : GCell()
  {
    cx_ = cx;
    cy_ = cy;
    dDx_ = dx_ = dx;
    dDy_ = dy_ = dy;
    densityScale_ = 1.0;
  }

  void GCell::setInstance(Instance *inst)
  {
    inst_ = inst;
    isMacro_ = inst->isMacro();
    double lx = static_cast<double>(inst->lx());
    double ly = static_cast<double>(inst->ly());
    double ux = static_cast<double>(inst->ux());
    double uy = static_cast<double>(inst->uy());
    cx_ = 0.5 * (lx + ux);
    cy_ = 0.5 * (ly + uy);
    dDx_ = dx_ = ux - lx;
    dDy_ = dy_ = uy - ly;
  }

  void GCell::addGPin(GPin *gPin)
  {
    gPins_.push_back(gPin);
  }

  void GCell::updatePins(bool useTheta)
  {
    if (useTheta)
    {
      for (GPin *pin : gPins_)
        pin->updateLocationWithTheta(this);
    }
    else
    {
      for (GPin *pin : gPins_)
        pin->updateLocation(this);
    }
  }

  ////////////////////////////////////////////////
  // GNet

  GNet::GNet()
      : net_(nullptr),
        lx_(0), ly_(0), ux_(0), uy_(0),
        waExpMinSumX_(0), waXExpMinSumX_(0),
        waExpMaxSumX_(0), waXExpMaxSumX_(0),
        waExpMinSumY_(0), waYExpMinSumY_(0),
        waExpMaxSumY_(0), waYExpMaxSumY_(0)
  {
  }

  GNet::GNet(Net *net)
      : GNet()
  {
    net_ = net;
  }

  void GNet::addGPin(GPin *gPin)
  {
    gPins_.push_back(gPin);
  }

  void GNet::updateBox()
  {
    lx_ = ly_ = std::numeric_limits<double>::max();
    ux_ = uy_ = std::numeric_limits<double>::lowest();

    for (auto &gPin : gPins_)
    {
      lx_ = std::min(gPin->cx(), lx_);
      ly_ = std::min(gPin->cy(), ly_);
      ux_ = std::max(gPin->cx(), ux_);
      uy_ = std::max(gPin->cy(), uy_);
    }
  }

  void GNet::clearWaVars()
  {
    waExpMinSumX_ = 0;
    waXExpMinSumX_ = 0;

    waExpMaxSumX_ = 0;
    waXExpMaxSumX_ = 0;

    waExpMinSumY_ = 0;
    waYExpMinSumY_ = 0;

    waExpMaxSumY_ = 0;
    waYExpMaxSumY_ = 0;
  }

  ////////////////////////////////////////////////
  // GPin

  GPin::GPin()
      : pin_(nullptr), gCell_(nullptr), gNet_(nullptr),
        offsetCx_(0), offsetCy_(0),
        cx_(0), cy_(0),
        maxExpSumX_(0), maxExpSumY_(0),
        minExpSumX_(0), minExpSumY_(0),
        hasMaxExpSumX_(0), hasMaxExpSumY_(0),
        hasMinExpSumX_(0), hasMinExpSumY_(0)
  {
  }

  GPin::GPin(Pin *pin)
      : GPin()
  {
    pin_ = pin;
    cx_ = static_cast<double>(pin->cx());
    cy_ = static_cast<double>(pin->cy());
    offsetCx_ = static_cast<double>(pin->offsetCx());
    offsetCy_ = static_cast<double>(pin->offsetCy());
  }

  void GPin::setCenterLocation(double cx, double cy)
  {
    cx_ = cx;
    cy_ = cy;
  }

  void GPin::clearWaVars()
  {
    hasMaxExpSumX_ = hasMaxExpSumY_ = false;
    hasMinExpSumX_ = hasMinExpSumY_ = false;
    maxExpSumX_ = maxExpSumY_ = 0;
    minExpSumX_ = minExpSumY_ = 0;
  }

  void GPin::updateLocation(const GCell *gCell)
  {
    cx_ = gCell->cx() + offsetCx_;
    cy_ = gCell->cy() + offsetCy_;
  }

  void GPin::updateLocationWithTheta(const GCell *gCell)
  {
    double cosTheta = std::cos(gCell->theta());
    double sinTheta = std::sin(gCell->theta());
    cx_ = gCell->cx() + offsetCx_ * cosTheta - offsetCy_ * sinTheta;
    cy_ = gCell->cy() + offsetCx_ * sinTheta + offsetCy_ * cosTheta;
  }

  ////////////////////////////////////////////////////////
  // Bin

  Bin::Bin()
      : lx_(0), ly_(0),
        ux_(0), uy_(0),
        nonPlaceArea_(0), instPlacedArea_(0),
        fillerArea_(0),
        density_(0),
        electroPhi_(0),
        electroForceX_(0), electroForceY_(0)
  {
  }

  Bin::Bin(double lx, double ly, double ux, double uy)
      : Bin()
  {
    lx_ = lx;
    ly_ = ly;
    ux_ = ux;
    uy_ = uy;
  }

  ////////////////////////////////////////////////
  // BinGrid

  BinGrid::BinGrid()
      : die_(nullptr),
        lx_(0), ly_(0), ux_(0), uy_(0),
        binCntX_(0), binCntY_(0),
        binSizeX_(0), binSizeY_(0),
        targetDensity_(0), overflowArea_(0),
        totalCellArea_(0)
  {
  }

  BinGrid::BinGrid(Die *die)
      : BinGrid()
  {
    setDie(die);
  }

  void BinGrid::setDie(Die *die)
  {
    die_ = die;
    lx_ = static_cast<double>(die->coreLx());
    ly_ = static_cast<double>(die->coreLy());
    ux_ = static_cast<double>(die->coreUx());
    uy_ = static_cast<double>(die->coreUy());

    placeInstCnt_ = 0;
    placeStdCellArea_ = 0;
    placeMacroArea_ = 0;
    fixedInstArea_ = 0;
    for (Instance *inst : die_->insts())
    {
      double area = static_cast<double>(inst->dx() * inst->dy());
      if (inst->isFixed())
        fixedInstArea_ += area;
      else
      {
        placeInstCnt_++;
        if (inst->isMacro())
          placeMacroArea_ += area;
        else
          placeStdCellArea_ += area;
      }
    }
  }

  void BinGrid::initBins()
  {
    int foundBinCnt = 2;
    if (placeInstCnt_ != 0)
    {
      // assume no fixed instance
      double totalBinArea = (double)(ux_ - lx_) * (uy_ - ly_);
      double avgPlaceInstArea = (placeMacroArea_ + placeStdCellArea_) / placeInstCnt_;
      double idealBinArea = avgPlaceInstArea / targetDensity_;
      int idealBinCnt = static_cast<int>(totalBinArea / idealBinArea);

      LOG_DEBUG("TargetDensity: {}", targetDensity_);
      LOG_DEBUG("AveragePlaceInstArea: {}", avgPlaceInstArea);
      LOG_DEBUG("IdealBinArea: {}", idealBinArea);
      LOG_DEBUG("IdealBinCnt: {}", idealBinCnt);
      LOG_DEBUG("TotalBinArea: {}", totalBinArea);

      // find binCnt: 2, 4, 8, 16, 32, 64, ...
      // s.t. binCnt^2 <= idealBinCnt <= (binCnt*2)^2.
      for (foundBinCnt = 2; foundBinCnt <= 1024; foundBinCnt *= 2)
      {
        if (foundBinCnt * foundBinCnt <= idealBinCnt && 4 * foundBinCnt * foundBinCnt > idealBinCnt)
          break;
      }
    }
    initBins(foundBinCnt, foundBinCnt);
  }

  void BinGrid::initBins(int cntX, int cntY)
  {
    binCntX_ = cntX;
    binCntY_ = cntY;
    binSizeX_ = static_cast<double>((ux_ - lx_) / binCntX_);
    binSizeY_ = static_cast<double>((uy_ - ly_) / binCntY_);
    LOG_DEBUG("BinCnt: {}, {}", binCntX_, binCntY_);
    LOG_DEBUG("BinSize: {}, {}", binSizeX_, binSizeY_);

    // initialize binStor_, bins_ vector
    binStor_.reserve(binCntX_ * binCntY_);
    bins_.reserve(binCntX_ * binCntY_);
    for (int idxY = 0; idxY < binCntY_; idxY++)
    {
      double ly = static_cast<double>(idxY * binSizeY_);
      double uy = ly + binSizeY_;
      for (int idxX = 0; idxX < binCntX_; idxX++)
      {
        double lx = static_cast<double>(idxX * binSizeX_);
        double ux = lx + binSizeX_;
        binStor_.emplace_back(lx, ly, ux, uy);
        bins_.push_back(&binStor_.back());
      }
    }

    // only initialized once
    updateBinsNonPlaceArea();

    // initialize fft structrue based on bins
    fft_.init(binCntX_, binCntY_, binSizeX_, binSizeY_);
  }

  void BinGrid::updateBinsNonPlaceArea()
  {
    for (auto &bin : bins_)
      bin->setNonPlaceArea(0);

    for (Instance *inst : die_->insts())
    {
      if (!inst->isFixed())
        continue;
      auto pairX = getMinMaxIdxX(static_cast<double>(inst->lx()),
                                 static_cast<double>(inst->ux()));
      auto pairY = getMinMaxIdxY(static_cast<double>(inst->ly()),
                                 static_cast<double>(inst->uy()));
      for (int j = pairY.first; j < pairY.second; j++)
      {
        for (int i = pairX.first; i < pairX.second; i++)
        {
          Bin *bin = bins_[j * binCntX_ + i];

          // Note that nonPlaceArea should have scale-down with
          // target density.
          // See MS-replace paper
          //
          double overlap = getOverlapArea(bin, inst);
          overlap *= targetDensity_;
          bin->addNonPlaceArea(overlap);
        }
      }
    }
  }

  // Core Part
  void BinGrid::updateBinsGCellDensityArea(bool useTheta)
  {
    // clear the Bin-area info
    for (auto &bin : bins_)
    {
      bin->setInstPlacedArea(0);
      bin->setFillerArea(0);
    }

    for (auto &cell : gCells_)
    {
      // The following function is critical runtime hotspot for global placer.
      if (cell->isMacro())
      {
        if (useTheta)
          addMacroBinAreaWithTheta(cell);
        else
          addMacroBinArea(cell);
      }
      else if (cell->isInstance())
      {
        addStdCellBinArea(cell);
      }
      else
      {
        addFillerBinArea(cell);
      }
    }

    overflowArea_ = 0;
    // update density and overflowArea
    // for nesterov use and FFT library
    for (auto &bin : bins_)
    {
      double binArea = bin->binArea();
      double cellArea = bin->instPlacedArea() + bin->fillerArea() + bin->nonPlaceArea();
      double density = cellArea / (binArea * targetDensity_);
      bin->setDensity(density);
      double binOverflowArea = bin->instPlacedArea() + bin->nonPlaceArea() - binArea * targetDensity_;
      overflowArea_ += std::max(binOverflowArea, 0.0);
    }
  }

  void BinGrid::addFillerBinArea(const GCell *gcell)
  {
    std::pair<int, int> pairX = getMinMaxIdxX(gcell->dLx(), gcell->dUx());
    std::pair<int, int> pairY = getMinMaxIdxY(gcell->dLy(), gcell->dUy());
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bins_[j * binCntX_ + i];
        double overlapArea = getOverlapDensityArea(bin, gcell);
        overlapArea *= gcell->densityScale();
        bin->addFillerArea(static_cast<double>(overlapArea));
      }
    }
  }

  void BinGrid::addStdCellBinArea(const GCell *gcell)
  {
    std::pair<int, int> pairX = getMinMaxIdxX(gcell->dLx(), gcell->dUx());
    std::pair<int, int> pairY = getMinMaxIdxY(gcell->dLy(), gcell->dUy());
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bins_[j * binCntX_ + i];
        double overlapArea = getOverlapDensityArea(bin, gcell);
        overlapArea *= gcell->densityScale();
        bin->addInstPlacedArea(static_cast<double>(overlapArea));
      }
    }
  }

  void BinGrid::addMacroBinArea(const GCell *gcell)
  {
    std::pair<int, int> pairX = getMinMaxIdxX(gcell->dLx(), gcell->dUx());
    std::pair<int, int> pairY = getMinMaxIdxY(gcell->dLy(), gcell->dUy());
    double scale = gcell->densityScale() * targetDensity_;
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bins_[j * binCntX_ + i];
        double overlapArea = getOverlapDensityArea(bin, gcell);
        overlapArea *= scale;
        bin->addInstPlacedArea(overlapArea);
      }
    }
  }

  void BinGrid::addMacroBinAreaWithTheta(const GCell *gcell)
  {
    std::pair<int, int> pairX = getMinMaxIdxX(gcell->dLx(), gcell->dUx());
    std::pair<int, int> pairY = getMinMaxIdxY(gcell->dLy(), gcell->dUy());
    double scale = gcell->densityScale() * bellShape1(gcell->theta()) * targetDensity_;
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bins_[j * binCntX_ + i];
        double overlapArea = getOverlapDensityArea(bin, gcell);
        overlapArea *= scale;
        bin->addInstPlacedArea(static_cast<double>(overlapArea));
      }
    }

    double rlx = gcell->cx() - 0.5f * gcell->dDy();
    double rux = gcell->cx() + 0.5f * gcell->dDx();
    double rly = gcell->cy() - 0.5f * gcell->dDy();
    double ruy = gcell->cy() + 0.5f * gcell->dDy();
    pairX = getMinMaxIdxX(rlx, rux);
    pairY = getMinMaxIdxY(rly, ruy);
    scale = gcell->densityScale() * bellShape2(gcell->theta()) * targetDensity_;
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bins_[j * binCntX_ + i];
        double overlapLx = std::max(rlx, bin->lx());
        double overlapUx = std::min(rux, bin->ux());
        double overlapLy = std::max(rly, bin->ly());
        double overlapUy = std::min(ruy, bin->uy());
        if (overlapLx < overlapUx && overlapLy < overlapUy)
        {
          double overlapArea = (overlapUx - overlapLx) * (overlapUy - overlapLy);
          overlapArea *= scale;
          bin->addInstPlacedArea(static_cast<double>(overlapArea));
        }
      }
    }
  }

  std::pair<int, int> BinGrid::getMinMaxIdxX(double lx1, double ux1)
  {
    int lowerIdx = static_cast<int>(std::floor((lx1 - lx_) / binSizeX_));
    int upperIdx = static_cast<int>(std::ceil((ux1 - lx_) / binSizeX_));
    return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntX_));
  }

  std::pair<int, int> BinGrid::getMinMaxIdxY(double ly1, double uy1)
  {
    int lowerIdx = static_cast<int>(std::floor((ly1 - ly_) / binSizeY_));
    int upperIdx = static_cast<int>(std::ceil((uy1 - ly_) / binSizeY_));
    return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntY_));
  }

  void BinGrid::addGCell(GCell *gCell)
  {
    gCells_.push_back(gCell);
    totalCellArea_ += gCell->dx() * gCell->dy();
  }

  void BinGrid::updateDensityForceBin()
  {
    // copy density to utilize FFT
    for (int y = 0; y < binCntY_; y++)
    {
      for (int x = 0; x < binCntX_; x++)
        fft_.updateDensity(x, y, bins_[y * binCntX_ + x]->density());
    }

    // do FFT
    fft_.doFFT();

    // update electroPhi and electroForce
    // update sumPhi_ for nesterov loop
    sumPhi_ = 0;
    for (int y = 0; y < binCntY_; y++)
    {
      for (int x = 0; x < binCntX_; x++)
      {
        auto eForcePair = fft_.getElectroForce(x, y);
        double electroPhi = fft_.getElectroPhi(x, y);

        Bin *bin = bins_[y * binCntX_ + x];
        bin->setElectroForceX(eForcePair.first);
        bin->setElectroForceY(eForcePair.second);
        bin->setElectroPhi(electroPhi);

        sumPhi_ += electroPhi * (bin->nonPlaceArea() + bin->instPlacedArea() + bin->fillerArea());
      }
    }
  }

  void BinGrid::updateGCellDensityScaleAndSize()
  {
    // update densitySize and densityScale in each gCell
    for (auto &gCell : gCells_)
    {
      double scaleX = 0, scaleY = 0;
      double densitySizeX = 0, densitySizeY = 0;
      if (gCell->dx() < SQRT2 * binSizeX_)
      {
        scaleX = static_cast<double>(gCell->dx()) / static_cast<double>(SQRT2 * binSizeX_);
        densitySizeX = static_cast<double>(SQRT2 * binSizeX_);
      }
      else
      {
        scaleX = 1.0;
        densitySizeX = gCell->dx();
      }

      if (gCell->dy() < SQRT2 * binSizeY_)
      {
        scaleY = static_cast<double>(gCell->dy()) / static_cast<double>(SQRT2 * binSizeY_);
        densitySizeY = static_cast<double>(SQRT2 * binSizeY_);
      }
      else
      {
        scaleY = 1.0;
        densitySizeY = gCell->dy();
      }

      // gCell->setSize(densitySizeX, densitySizeY);
      gCell->setDensitySize(densitySizeX, densitySizeY);
      gCell->setDensityScale(scaleX * scaleY);
    }
  }

  ////////////////////////////////////////////////
  // NesterovBaseVars
  NesterovBaseVars::NesterovBaseVars()
      : targetDensity(1.0),
        minAvgCut(0.1), maxAvgCut(0.9),
        binCntX(0), binCntY(0),
        minWireLengthForceBar(-300),
        isSetBinCntX(0), isSetBinCntY(0)
  {
  }

  void NesterovBaseVars::reset()
  {
    targetDensity = 1.0;
    minAvgCut = 0.1;
    maxAvgCut = 0.9;
    isSetBinCntX = isSetBinCntY = 0;
    binCntX = binCntY = 0;
    minWireLengthForceBar = -300;
  }

  ////////////////////////////////////////////////
  // NesterovBase

  NesterovBase::NesterovBase()
      : pb_(nullptr), sumPhi_(), nbVars_()
  {
  }

  NesterovBase::NesterovBase(NesterovBaseVars nbVars, std::shared_ptr<PlacerBase> pb)
      : NesterovBase()
  {
    nbVars_ = nbVars;
    pb_ = pb;
    init();
  }

  void NesterovBase::init()
  {
    LOG_TRACE("start NesterovBase::init");
    std::unordered_map<Instance *, GCell *> gCellMap;
    std::unordered_map<Pin *, GPin *> gPinMap;
    std::unordered_map<Net *, GNet *> gNetMap;
    std::unordered_map<Die *, BinGrid *> binGridMap;

    // rearrange instances
    // |Fixed Cell|Movable Macros|Moveble StdCells|
    std::vector<Instance *> &insts = pb_->insts();
    auto fixedIt = std::partition(insts.begin(), insts.end(), [](const Instance *inst)
                                  { return inst->isFixed(); });
    auto macroIt = std::partition(fixedIt, insts.end(), [](const Instance *inst)
                                  { return inst->isMacro(); });

    // gCellStor init
    gCellStor_.reserve(pb_->insts().size());
    for (Instance *inst : insts)
    {
      if (inst->isFixed())
        continue;
      GCell myGCell(inst);
      gCellStor_.push_back(myGCell);
    }

    // create filler
    std::unordered_map<Die*, std::pair<int, int>> dieFillerIdxMap;
    for(Die* die : pb_->dies())
    {
      int idx1 = static_cast<int>(gCellStor_.size());
      initFillerGCells(die);
      int idx2 = static_cast<int>(gCellStor_.size());
      dieFillerIdxMap.emplace(die, std::make_pair(idx1, idx2));
    }

    // generate gCell ptr
    gCells_.reserve(gCellStor_.size());
    gCellMap.reserve(gCellStor_.size());
    for (auto &gCell : gCellStor_)
    {
      gCells_.push_back(&gCell);
      if (gCell.isInstance())
        gCellMap.emplace(gCell.instance(), &gCell);
    }

    // create binGrid and generate binGrid ptr
    binGridStor_.reserve(pb_->dies().size());
    binGrids_.reserve(pb_->dies().size());
    binGridMap.reserve(pb_->dies().size());
    for (Die *die : pb_->dies())
    {
      // each die should own one bin grid
      binGridStor_.emplace_back(die);
      BinGrid* bg = &binGridStor_.back();
      binGrids_.push_back(bg);
      binGridMap.emplace(die, bg);

      // send param into binGrid structure
      bg->setTargetDensity(nbVars_.targetDensity);
      if (nbVars_.isSetBinCntX && nbVars_.isSetBinCntY)
        bg->initBins(nbVars_.binCntX, nbVars_.binCntY);
      else
        bg->initBins();

      // assign gCell to binGrid
      for (Instance *inst : die->insts())
      {
        if(inst->isFixed())
          continue;
        auto it = gCellMap.find(inst);
        GCell *gCell = (it == gCellMap.end() ? nullptr : it->second);
        bg->addGCell(gCell);
        gCell->setBinGrid(bg);
      }

      // assign filler to binGrid
      auto pair = dieFillerIdxMap.at(die);
      for(int i = pair.first; i < pair.second; i++)
      {
        bg->addGCell(gCells_[i]);
        gCells_[i]->setBinGrid(bg);
      }

      bg->updateGCellDensityScaleAndSize();
    }

    // gPinStor init and generate gPin ptr
    gPinStor_.reserve(pb_->pins().size());
    gPins_.reserve(pb_->pins().size());
    gPinMap.reserve(pb_->pins().size());
    for (Pin* pin : pb_->pins())
    {
      gPinStor_.emplace_back(pin);
      gPins_.push_back(&gPinStor_.back());
      gPinMap.emplace(pin, &gPinStor_.back());
    }

    // gNetStor init
    gNetStor_.reserve(pb_->nets().size());
    gNets_.reserve(pb_->nets().size());
    gNetMap.reserve(pb_->nets().size());
    for (Net* net : pb_->nets())
    {
      gNetStor_.emplace_back(net);
      gNets_.push_back(&gNetStor_.back());
      gNetMap.emplace(net, &gNetStor_.back());
    }

    // gCellStor_'s pins_ fill
    for (GCell* gCell : gCells_)
    {
      if (gCell->isInstance())
      {
        for (auto &pin : gCell->instance()->pins())
        {
          auto it = gPinMap.find(pin);
          gCell->addGPin(it == gPinMap.end() ? nullptr : it->second);
        }
      }
    }

    // gPinStor_' GNet and GCell fill
    for (GPin* gPin : gPins_)
    {
      auto iit = gCellMap.find(gPin->pin()->instance());
      auto nit = gNetMap.find(gPin->pin()->net());
      gPin->setGCell(iit == gCellMap.end() ? nullptr : iit->second);
      gPin->setGNet(nit == gNetMap.end() ? nullptr : nit->second);
    }

    // gNetStor_'s GPin fill
    for (GNet* gNet : gNets_)
    {
      for (Pin* pin : gNet->net()->pins())
      {
        // NOTE: What if a pin is on a fixed instances?
        auto it = gPinMap.find(pin);
        gNet->addGPin(it == gPinMap.end() ? nullptr : it->second);
      }
    }

    LOG_DEBUG("FillerInit: NumGCells: {}", gCells_.size());
    LOG_DEBUG("FillerInit: NumGNets: {}", gNets_.size());
    LOG_DEBUG("FillerInit: NumGPins: {}", gPins_.size());
    LOG_TRACE("finish NesterovBase::init");
  }

  // virtual filler GCells
  void NesterovBase::initFillerGCells(const Die *die)
  {
    double dxSum = 0.0, dySum = 0.0;
    int numMovableCell = 0;
    double fixedInstanceArea = 0;
    double placeStdCellArea = 0;
    double placeMacroArea = 0;
    for (const Instance *inst : die->insts())
    {
      double area = static_cast<double>(inst->dx()) * inst->dy();
      if(inst->isFixed())
        fixedInstanceArea += area;
      else if(inst->isMacro())
        placeMacroArea += area;
      else
        placeStdCellArea += area;
      
      if(!inst->isFixed())
      {
        dxSum += inst->dx();
        dySum += inst->dy();
        numMovableCell++;
      }
    }
    // if die doesn't have any moveable instance, we will not create any filler.
    if(numMovableCell == 0)
      return;

    double avgDx = dxSum / numMovableCell;
    double avgDy = dySum / numMovableCell;
    double coreArea = static_cast<double>(die->coreDx()) * die->coreDy();
    // nonPlaceInstsArea should not have targetDensity downscaling!!!
    double whiteSpaceArea = coreArea - fixedInstanceArea;
    // TODO density screening
    double movableArea = whiteSpaceArea * nbVars_.targetDensity;
    double totalFillerArea = movableArea - placeStdCellArea - placeMacroArea * nbVars_.targetDensity;
    if (totalFillerArea < 0)
    {
      LOG_ERROR("Filler area is negative!!\n"
                "\tPlease put higher target density or\n"
                "\tRe-floorplan to have enough coreArea");
    }

    int fillerCnt = static_cast<int>(totalFillerArea / (avgDx * avgDy));
    LOG_DEBUG("FillerInit: CoreArea: {}", coreArea);
    LOG_DEBUG("FillerInit: WhiteSpaceArea: {}", whiteSpaceArea);
    LOG_DEBUG("FillerInit: MovableArea: {}", movableArea);
    LOG_DEBUG("FillerInit: TotalFillerArea: {}", totalFillerArea);
    LOG_DEBUG("FillerInit: NumFillerCells: {}", fillerCnt);
    LOG_DEBUG("FillerInit: FillerCellArea: {}", avgDx * avgDy);
    LOG_DEBUG("FillerInit: FillerCellSize: {}, {}", avgDx, avgDy);

    //
    // mt19937 supports huge range of random values.
    // rand()'s RAND_MAX is only 32767.
    //
    std::mt19937 randVal(0);
    for (int i = 0; i < fillerCnt; i++)
    {
      // instability problem between g++ and clang++!
      auto randX = randVal();
      auto randY = randVal();

      // place filler cells on random coordi and
      // set size as avgDx and avgDy
      double cx = static_cast<double>(randX % die->coreDx() + die->coreLx());
      double cy = static_cast<double>(randY % die->coreDy() + die->coreLy());
      GCell myGCell(cx, cy, avgDx, avgDy);
      gCellStor_.emplace_back(cx, cy, avgDx, avgDy);
    }
  }

  // gcell update
  void NesterovBase::updateGCellCenterLocation(const std::vector<Point> &coordis)
  {
    for(int i = 0; i < coordis.size(); i++)
    {
      gCells_[i]->setCenterLocation(coordis[i].x, coordis[i].y);
      gCells_[i]->updatePins(false);
    }
    for (BinGrid *bg : binGrids_)
    {
      bg->updateBinsGCellDensityArea(false);
    }
  }

  void NesterovBase::updateGCellCenterLocationWithTheta(const std::vector<Point> &coordis, const std::vector<double> &macrosTheta)
  {
    // the first macroCount gCells ara macros;
    int macroCount = static_cast<int>(macrosTheta.size());
    for(int i = 0; i < macroCount; i++)
    {
      gCells_[i]->setTheta(macrosTheta[i]);
      gCells_[i]->setCenterLocation(coordis[i].x, coordis[i].y);
      gCells_[i]->updatePins(true);
    }
    for(int i = macroCount; i < coordis.size(); i++)
    {
      gCells_[i]->setCenterLocation(coordis[i].x, coordis[i].y);
      gCells_[i]->updatePins(false);
    }
    for (BinGrid *bg : binGrids_)
    {
      bg->updateBinsGCellDensityArea(true);
    }
  }

  double NesterovBase::overflowArea() const
  {
    double area = 0;
    for (BinGrid *bg : binGrids_)
      area += bg->overflowArea();
    return area;
  }

  double NesterovBase::sumPhi() const
  {
    return sumPhi_;
  }

  double NesterovBase::targetDensity() const
  {
    return nbVars_.targetDensity;
  }

  void NesterovBase::updateDensityCoordiLayoutInside(GCell *gcell, bool useTheta)
  {
    double cx = getDensityCoordiLayoutInsideX(gcell, gcell->cx(), useTheta);
    double cy = getDensityCoordiLayoutInsideY(gcell, gcell->cy(), useTheta);
    gcell->setCenterLocation(cx, cy);
    gcell->updatePins(false);
  }

  double NesterovBase::getDensityCoordiLayoutInsideX(const GCell *gcell, double newCx, bool useTheta)
  {
    double adjVal = newCx;
    const BinGrid *bg = gcell->binGrid();
    adjVal = std::max(adjVal, bg->lx() + 0.5f * gcell->dDx());
    adjVal = std::min(adjVal, bg->ux() - 0.5f * gcell->dDx());
    if (useTheta)
    {
      adjVal = std::max(adjVal, bg->lx() + 0.5f * gcell->dDy());
      adjVal = std::min(adjVal, bg->ux() - 0.5f * gcell->dDy());
    }
    return adjVal;
  }

  double NesterovBase::getDensityCoordiLayoutInsideY(const GCell *gcell, double newCy, bool useTheta)
  {
    double adjVal = newCy;
    const BinGrid *bg = gcell->binGrid();
    adjVal = std::max(adjVal, bg->ly() + 0.5f * gcell->dDy());
    adjVal = std::min(adjVal, bg->uy() - 0.5f * gcell->dDy());
    if (useTheta)
    {
      adjVal = std::max(adjVal, bg->ly() + 0.5f * gcell->dDx());
      adjVal = std::min(adjVal, bg->uy() - 0.5f * gcell->dDx());
    }
    return adjVal;
  }

  //
  // WA force cals - wlCoeffX / wlCoeffY
  //
  // * Note that wlCoeffX and wlCoeffY is 1/gamma
  // in ePlace paper.
  void NesterovBase::updateWireLengthForceWA(double wlCoeffX, double wlCoeffY)
  {
    // clear all WA variables.
    for (auto &gNet : gNets_)
      gNet->clearWaVars();
    for (auto &gPin : gPins_)
      gPin->clearWaVars();

    for (auto &gNet : gNets_)
    {
      gNet->updateBox();

      for (auto &gPin : gNet->gPins())
      {
        double expMinX = (gNet->lx() - gPin->cx()) * wlCoeffX;
        double expMaxX = (gPin->cx() - gNet->ux()) * wlCoeffX;
        double expMinY = (gNet->ly() - gPin->cy()) * wlCoeffY;
        double expMaxY = (gPin->cy() - gNet->uy()) * wlCoeffY;

        // min x
        if (expMinX > nbVars_.minWireLengthForceBar)
        {
          gPin->setMinExpSumX(fastExp(expMinX));
          gNet->addWaExpMinSumX(gPin->minExpSumX());
          gNet->addWaXExpMinSumX(gPin->cx() * gPin->minExpSumX());
        }

        // max x
        if (expMaxX > nbVars_.minWireLengthForceBar)
        {
          gPin->setMaxExpSumX(fastExp(expMaxX));
          gNet->addWaExpMaxSumX(gPin->maxExpSumX());
          gNet->addWaXExpMaxSumX(gPin->cx() * gPin->maxExpSumX());
        }

        // min y
        if (expMinY > nbVars_.minWireLengthForceBar)
        {
          gPin->setMinExpSumY(fastExp(expMinY));
          gNet->addWaExpMinSumY(gPin->minExpSumY());
          gNet->addWaYExpMinSumY(gPin->cy() * gPin->minExpSumY());
        }

        // max y
        if (expMaxY > nbVars_.minWireLengthForceBar)
        {
          gPin->setMaxExpSumY(fastExp(expMaxY));
          gNet->addWaExpMaxSumY(gPin->maxExpSumY());
          gNet->addWaYExpMaxSumY(gPin->cy() * gPin->maxExpSumY());
        }
      }
    }
  }

  // get x,y WA Gradient values from GPin
  // Please check the JingWei's Ph.D. thesis full paper,
  // Equation (4.13)
  //
  // You can't understand the following function
  // unless you read the (4.13) formula
  Point NesterovBase::getWireLengthGradientPinWA(const GPin *gPin, double wlCoeffX, double wlCoeffY)
  {
    double gradientMinX = 0, gradientMinY = 0;
    double gradientMaxX = 0, gradientMaxY = 0;

    // min x
    if (gPin->hasMinExpSumX())
    {
      // from Net.
      double waExpMinSumX = gPin->gNet()->waExpMinSumX();
      double waXExpMinSumX = gPin->gNet()->waXExpMinSumX();

      gradientMinX =
          (waExpMinSumX * (gPin->minExpSumX() * (1.0 - wlCoeffX * gPin->cx())) + wlCoeffX * gPin->minExpSumX() * waXExpMinSumX) / (waExpMinSumX * waExpMinSumX);
    }

    // max x
    if (gPin->hasMaxExpSumX())
    {

      double waExpMaxSumX = gPin->gNet()->waExpMaxSumX();
      double waXExpMaxSumX = gPin->gNet()->waXExpMaxSumX();

      gradientMaxX =
          (waExpMaxSumX * (gPin->maxExpSumX() * (1.0 + wlCoeffX * gPin->cx())) - wlCoeffX * gPin->maxExpSumX() * waXExpMaxSumX) / (waExpMaxSumX * waExpMaxSumX);
    }

    // min y
    if (gPin->hasMinExpSumY())
    {

      double waExpMinSumY = gPin->gNet()->waExpMinSumY();
      double waYExpMinSumY = gPin->gNet()->waYExpMinSumY();

      gradientMinY =
          (waExpMinSumY * (gPin->minExpSumY() * (1.0 - wlCoeffY * gPin->cy())) + wlCoeffY * gPin->minExpSumY() * waYExpMinSumY) / (waExpMinSumY * waExpMinSumY);
    }

    // max y
    if (gPin->hasMaxExpSumY())
    {

      double waExpMaxSumY = gPin->gNet()->waExpMaxSumY();
      double waYExpMaxSumY = gPin->gNet()->waYExpMaxSumY();

      gradientMaxY =
          (waExpMaxSumY * (gPin->maxExpSumY() * (1.0 + wlCoeffY * gPin->cy())) - wlCoeffY * gPin->maxExpSumY() * waYExpMaxSumY) / (waExpMaxSumY * waExpMaxSumY);
    }

    return Point(gradientMinX - gradientMaxX,
                 gradientMinY - gradientMaxY);
  }

  // get x,y WA Gradient values with given GCell
  Point NesterovBase::getWireLengthGradientWA(const GCell *gCell, double wlCoeffX, double wlCoeffY)
  {
    Point gradientPair;
    for (auto &gPin : gCell->gPins())
    {
      auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);
      gradientPair.x += tmpPair.x;
      gradientPair.y += tmpPair.y;
    }
    return gradientPair;
  }

  Point NesterovBase::getWireLengthGradientWAWithTheta(const GCell *gCell, double wlCoeffX, double wlCoeffY, double &gradTheta)
  {
    Point gradientPair;
    gradTheta = 0.0;
    double cosTheta = std::cos(gCell->theta());
    double sinTheta = std::sin(gCell->theta());
    for (auto &gPin : gCell->gPins())
    {
      auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);
      gradientPair.x += tmpPair.x;
      gradientPair.y += tmpPair.y;
      gradTheta += tmpPair.x * (-gPin->offsetCx() * sinTheta - gPin->offsetCy() * cosTheta);
      gradTheta += tmpPair.y * (gPin->offsetCx() * cosTheta - gPin->offsetCy() * sinTheta);
    }
    return gradientPair;
  }

  // get GCells' electroForcePair
  // i.e. get DensityGradient with given GCell
  Point NesterovBase::getDensityGradient(const GCell *gCell)
  {
    BinGrid *bg = gCell->binGrid();
    auto pairX = bg->getMinMaxIdxX(gCell->dLx(), gCell->dUx());
    auto pairY = bg->getMinMaxIdxY(gCell->dLy(), gCell->dUy());
    Point electroForce;
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bg->bins()[j * bg->binCntX() + i];
        double overlapArea = getOverlapDensityArea(bin, gCell) * gCell->densityScale();

        electroForce.x += overlapArea * bin->electroForceX();
        electroForce.y += overlapArea * bin->electroForceY();
      }
    }
    return electroForce;
  }

  Point NesterovBase::getDensityGradientWithTheta(const GCell *gCell, double &gradTheta)
  {
    BinGrid *bg = gCell->binGrid();
    Point electroForce;
    gradTheta = 0.0;

    auto pairX = bg->getMinMaxIdxX(gCell->dLx(), gCell->dUx());
    auto pairY = bg->getMinMaxIdxY(gCell->dLy(), gCell->dUy());
    double bell = bellShape1(gCell->theta());
    double derive = bellShapeDerive1(gCell->theta());
    for (int i = pairX.first; i < pairX.second; i++)
    {
      for (int j = pairY.first; j < pairY.second; j++)
      {
        Bin *bin = bg->bins()[j * bg->binCntX() + i];
        double overlapArea = getOverlapDensityArea(bin, gCell) * gCell->densityScale();
        electroForce.x += bell * overlapArea * bin->electroForceX();
        electroForce.y += bell * overlapArea * bin->electroForceY();
        gradTheta -= derive * overlapArea * bin->electroPhi();
      }
    }

    double rlx = gCell->cx() - 0.5f * gCell->dDy();
    double rux = gCell->cx() + 0.5f * gCell->dDy();
    double rly = gCell->cy() - 0.5f * gCell->dDx();
    double ruy = gCell->cy() + 0.5f * gCell->dDx();
    pairX = bg->getMinMaxIdxX(rlx, rux);
    pairY = bg->getMinMaxIdxY(rly, ruy);
    bell = bellShape2(gCell->theta());
    derive = bellShapeDerive2(gCell->theta());
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bg->bins()[j * bg->binCntX() + i];
        double overlapLx = std::max(rlx, bin->lx());
        double overlapLy = std::max(rly, bin->ly());
        double overlapUx = std::min(rux, bin->ux());
        double overlapUy = std::min(ruy, bin->uy());
        if (overlapLx < overlapUx && overlapLy < overlapUy)
        {
          double overlapArea = (overlapUx - overlapLx) * (overlapUy - overlapLy) * gCell->densityScale();
          electroForce.x += bell * overlapArea * bin->electroForceX();
          electroForce.y += bell * overlapArea * bin->electroForceY();
          gradTheta -= derive * overlapArea * bin->electroPhi();
        }
      }
    }
    return electroForce;
  }

  Point NesterovBase::getLocalDensityGradient(const GCell *gCell, double alpha, double beta, double &cellDelta)
  {
    BinGrid *bg = gCell->binGrid();
    auto pairX = bg->getMinMaxIdxX(gCell->dLx(), gCell->dUx());
    auto pairY = bg->getMinMaxIdxY(gCell->dLy(), gCell->dUy());
    Point grad;

    double binArea = bg->binSizeX() * bg->binSizeY();
    double invTotalCellArea = static_cast<double>(1.0 / bg->totalCellArea());

    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bg->bins()[j * bg->binCntX() + i];
        double overlapLx = std::max(bin->lx(), gCell->dLx());
        double overlapLy = std::max(bin->ly(), gCell->dLy());
        double overlapUx = std::min(bin->ux(), gCell->dUx());
        double overlapUy = std::min(bin->uy(), gCell->dUy());
        if (overlapLx < overlapUx && overlapLy < overlapUy)
        {
          double binCellArea = bin->density() * binArea;
          if (binCellArea > binArea)
            cellDelta += beta * (binCellArea - binArea) * invTotalCellArea;

          double dAdx = 0.0f, dAdy = 0.0f;
          if (bin->lx() < gCell->dLx() && bin->ux() < gCell->dUx()) // ov_x = bin_ux - cell_lx
            dAdx = -(overlapUy - overlapLy);
          else if (bin->lx() > gCell->dLx() && bin->ux() > gCell->dUx()) // ov_x = cell_ux - bin_lx
            dAdx = (overlapUy - overlapLy);
          if (bin->ly() < gCell->dLy() && bin->uy() < gCell->dUy())
            dAdy = -(overlapUx - overlapLx);
          else if (bin->ly() > gCell->dLy() && bin->uy() > gCell->dUy())
            dAdy = (overlapUx - overlapLx);
            
          double shareArea = gCell->densityScale() * (overlapUx - overlapLx) * (overlapUy - overlapLy);
          double commonDiv = alpha / binArea;
          double nu = fastExp(commonDiv * (binCellArea - binArea));
          double commonVal = commonDiv * nu * binCellArea * gCell->densityScale() * bin->electroPhi();
          double commonMul = nu * shareArea;
          grad.x += dAdx * commonVal;
          grad.x += commonMul * bin->electroForceX();
          grad.y += dAdy * commonVal;
          grad.y += commonMul * bin->electroForceY();
        }
      }
    }
    return grad;
  }

  Point NesterovBase::getLocalDensityGradientWithTheta(const GCell *gCell, double alpha, double beta, double &cellDelta, double &gradTheta)
  {
    BinGrid *bg = gCell->binGrid();
    auto pairX = bg->getMinMaxIdxX(gCell->dLx(), gCell->dUx());
    auto pairY = bg->getMinMaxIdxY(gCell->dLy(), gCell->dUy());
    Point grad;
    gradTheta = 0.0f;

    double binArea = bg->binSizeX() * bg->binSizeY();
    double invTotalCellArea = static_cast<double>(1.0 / bg->totalCellArea());
    double bell = bellShape1(gCell->theta());
    double derive = bellShapeDerive1(gCell->theta());
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bg->bins()[j * bg->binCntX() + i];
        double overlapLx = std::max(bin->lx(), gCell->dLx());
        double overlapLy = std::max(bin->ly(), gCell->dLy());
        double overlapUx = std::min(bin->ux(), gCell->dUx());
        double overlapUy = std::min(bin->uy(), gCell->dUy());
        if (overlapLx < overlapUx && overlapLy < overlapUy)
        {
          double binCellArea = bin->density() * binArea;
          if (binCellArea > binArea)
            cellDelta += beta * (binCellArea - binArea) * invTotalCellArea;

          double dAdx = 0.0f, dAdy = 0.0f;
          if (bin->lx() < gCell->dLx() && bin->ux() < gCell->dUx())
            dAdx = -(overlapUy - overlapLy);
          else if (bin->lx() > gCell->dLx() && bin->ux() > gCell->dUx())
            dAdx = (overlapUy - overlapLy);
          if (bin->ly() < gCell->dLy() && bin->uy() < gCell->dUy())
            dAdy = -(overlapUx - overlapLx);
          else if (bin->ly() > gCell->dLy() && bin->uy() > gCell->dUy())
            dAdy = (overlapUx - overlapLx);

          double shareArea = gCell->densityScale() * (overlapUx - overlapLx) * (overlapUy - overlapLy);
          double commonDiv = alpha / binArea;
          double nu = fastExp(commonDiv * (binCellArea - binArea));
          double commonVal = commonDiv * nu * binCellArea * gCell->densityScale() * bin->electroPhi();
          double commonMul = nu * shareArea;
          grad.x += bell * dAdx * commonVal;
          grad.x += bell * commonMul * bin->electroForceX();
          grad.y += bell * dAdy * commonVal;
          grad.y += bell * commonMul * bin->electroForceY();
          gradTheta -= derive * nu * shareArea * bin->electroPhi();
        }
      }
    }

    double rlx = gCell->cx() - 0.5f * gCell->dDy();
    double rux = gCell->cx() + 0.5f * gCell->dDy();
    double rly = gCell->cy() - 0.5f * gCell->dDx();
    double ruy = gCell->cy() + 0.5f * gCell->dDx();
    pairX = bg->getMinMaxIdxX(rlx, rux);
    pairY = bg->getMinMaxIdxY(rly, ruy);
    bell = bellShape2(gCell->theta());
    derive = bellShapeDerive2(gCell->theta());
    for (int j = pairY.first; j < pairY.second; j++)
    {
      for (int i = pairX.first; i < pairX.second; i++)
      {
        Bin *bin = bg->bins()[j * bg->binCntX() + i];
        double overlapLx = std::max(bin->lx(), rlx);
        double overlapLy = std::max(bin->ly(), rly);
        double overlapUx = std::min(bin->ux(), rux);
        double overlapUy = std::min(bin->uy(), ruy);
        if (overlapLx < overlapUx && overlapLy < overlapUy)
        {
          double binCellArea = bin->density() * binArea;
          if (binCellArea > binArea)
            cellDelta += beta * (binCellArea - binArea) * invTotalCellArea;

          double dAdx = 0.0f, dAdy = 0.0f;
          if (bin->lx() < rlx && bin->ux() < rux)
            dAdx = -(overlapUy - overlapLy);
          else if (bin->lx() > rlx && bin->ux() > rux)
            dAdx = (overlapUy - overlapLy);
          if (bin->ly() < rly && bin->uy() < ruy)
            dAdy = -(overlapUx - overlapLx);
          else if (bin->ly() > rly && bin->uy() > ruy)
            dAdy = (overlapUx - overlapLx);

          double shareArea = gCell->densityScale() * (overlapUx - overlapLx) * (overlapUy - overlapLy);
          double commonDiv = alpha / binArea;
          double nu = fastExp(commonDiv * (binCellArea - binArea));
          double commonVal = commonDiv * nu * binCellArea * gCell->densityScale() * bin->electroPhi();
          double commonMul = nu * shareArea;
          grad.x += bell * dAdx * commonVal;
          grad.x += bell * commonMul * bin->electroForceX();
          grad.y += bell * dAdy * commonVal;
          grad.y += bell * commonMul * bin->electroForceY();
          gradTheta -= derive * nu * shareArea * bin->electroPhi();
        }
      }
    }
    return grad;
  }

  Point NesterovBase::getWireLengthPreconditioner(const GCell *gCell)
  {
    double sizeVal = static_cast<double>(gCell->gPins().size());
    return Point(sizeVal, sizeVal);
  }

  Point NesterovBase::getDensityPreconditioner(const GCell *gCell)
  {
    // double areaVal = gCell->dDx() * gCell->dDy();
    double areaVal = gCell->dx() * gCell->dy();
    return Point(areaVal, areaVal);
  }

  Point NesterovBase::getLocalDensityPreconditioner(const GCell *gCell)
  {
    // double areaVal = gCell->dDx() * gCell->dDy();
    double areaVal = gCell->dx() * gCell->dy();
    return Point(areaVal, areaVal);
  }

  double NesterovBase::getWireLengthPreconditionerTheta(const GCell *gCell)
  {
    double precondi = 0;
    for (const GPin *gPin : gCell->gPins())
      precondi += gPin->offsetCx() * gPin->offsetCx() + gPin->offsetCy() * gPin->offsetCy();
    precondi *= 0.5;
    return precondi;
  }

  double NesterovBase::getDensityPreconditionerTheta(const GCell *gCell)
  {
    // return gCell->dDx() * gCell->dDy();
    return gCell->dx() * gCell->dy();
  }

  double NesterovBase::getLocalDensityPreconditionerTheta(const GCell *gCell)
  {
    // return gCell->dDx() * gCell->dDy();
    return gCell->dx() * gCell->dy();
  }

  // Density force cals
  void NesterovBase::updateDensityForceBin()
  {
    sumPhi_ = 0.0f;
    for (BinGrid *bg : binGrids_)
    {
      bg->updateDensityForceBin();
      sumPhi_ += bg->sumPhi();
    }
  }

  void NesterovBase::updateNetsBox()
  {
    for (auto &gNet : gNets_)
      gNet->updateBox();
  }

  double NesterovBase::getHpwl() const
  {
    double hpwl = 0;
    for (auto &gNet : gNets_)
      hpwl += gNet->hpwl();
    return hpwl;
  }

  double NesterovBase::getOverflow() const
  {
    double overflowArea = 0;
    double placeStdCellArea = 0;
    double placeMacroArea = 0;
    for (BinGrid *bg : binGrids_)
    {
      overflowArea += bg->overflowArea();
      placeStdCellArea += bg->placeStdCellArea();
      placeMacroArea += bg->placeMacroArea();
    }
    return overflowArea / (placeStdCellArea + placeMacroArea * nbVars_.targetDensity);
  }

  // int64_t is recommended, but double is 2x fast
  static double getOverlapDensityArea(const Bin *bin, const GCell *cell)
  {
    double lx = std::max(bin->lx(), cell->dLx()),
           ly = std::max(bin->ly(), cell->dLy()),
           ux = std::min(bin->ux(), cell->dUx()),
           uy = std::min(bin->uy(), cell->dUy());

    if (lx < ux && ly < uy)
      return (ux - lx) * (uy - ly);
    else
      return 0.0;
  }

  static double getOverlapArea(const Bin *bin, const Instance *inst)
  {
    double lx = std::max(bin->lx(), static_cast<double>(inst->lx())),
           ly = std::max(bin->ly(), static_cast<double>(inst->ly())),
           ux = std::min(bin->ux(), static_cast<double>(inst->ux())),
           uy = std::min(bin->uy(), static_cast<double>(inst->uy()));

    if (lx < ux && ly < uy)
      return (ux - lx) * (uy - ly);
    else
      return 0.0;
  }

  //
  // https://codingforspeed.com/using-faster-exponential-approximation/
  static double fastExp(double a)
  {
    a = 1.0 + a / 1024.0;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    return a;
  }

  static double bellShape(double theta_i, double theta)
  {
    double delta_abs = (1.0 / (2.0 * PI)) * std::abs(theta_i - theta);
    if (delta_abs <= 0.125)
      return 1.0 - 32.0 * delta_abs * delta_abs;
    else if (delta_abs <= 0.25)
      return 32.0 * (delta_abs - 0.25) * (delta_abs - 0.25);
    else
      return 0.0;
  }

  static double bellShape1(double theta_i)
  {
    return bellShape(theta_i, 0.0) + bellShape(theta_i, PI) + bellShape(theta_i, 2.0 * PI);
  }

  static double bellShape2(double theta_i)
  {
    return bellShape(theta_i, 0.5 * PI) + bellShape(theta_i, 1.5 * PI);
  }

  static double bellShapeDerive(double theta_i, double theta)
  {
    double delta = (1.0 / (2.0 * PI)) * (theta_i - theta);
    double delta_abs = std::abs(delta);
    if (delta_abs <= 0.125)
      return (-32.0 / PI) * delta;
    else if (delta_abs <= 0.25)
      return (32.0 / PI) * delta - std::copysign(8.0 / PI, delta);
    else
      return 0.0;
  }

  static double bellShapeDerive1(double theta_i)
  {
    return bellShapeDerive(theta_i, 0.0) + bellShapeDerive(theta_i, PI) + bellShapeDerive(theta_i, 2.0 * PI);
  }

  static double bellShapeDerive2(double theta_i)
  {
    return bellShapeDerive(theta_i, 0.5 * PI) + bellShapeDerive(theta_i, 1.5 * PI);
  }

}
