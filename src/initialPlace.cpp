#include "initialPlace.h"
#include "placerBase.h"
#include "log.h"
#include "plot.h"

#include <Eigen/IterativeLinearSolvers>

namespace replace
{
  using namespace std;

  using Eigen::BiCGSTAB;
  using Eigen::IdentityPreconditioner;

  typedef Eigen::Triplet<float> T;

  InitialPlaceVars::InitialPlaceVars()
      : maxIter(20),
        minDiffLength(1500),
        maxSolverIter(100),
        maxFanout(200),
        netWeightScale(800.0),
        incrementalPlaceMode(false)
  {
  }

  void InitialPlaceVars::reset()
  {
    maxIter = 20;
    minDiffLength = 1500;
    maxSolverIter = 100;
    maxFanout = 200;
    netWeightScale = 800.0;
    incrementalPlaceMode = false;
  }

  InitialPlace::InitialPlace()
      : ipVars_(), pb_(nullptr), placeInstCnt_(0)
  {
  }

  InitialPlace::InitialPlace(InitialPlaceVars ipVars, std::shared_ptr<PlacerBase> pb)
      : ipVars_(ipVars), pb_(pb)
  {
    placeInstCnt_ = 0;
    for(Instance* inst : pb->insts())
    {
      if(!inst->isFixed())
        placeInstCnt_++;
    }
  }

  void InitialPlace::doBicgstabPlace()
  {
    LOG_TRACE("Initial Placement Begin");

    float errorX = 0.0f, errorY = 0.0f;

    // normally, initial place will place all cells in the centers.
    if (!ipVars_.incrementalPlaceMode)
    {
      placeInstsCenter();
    }

    Plot::plot(pb_.get(), "./plot/cell", "before_ip");

    // set ExtId for idx reference // easy recovery
    setPlaceInstExtId();
    for (int i = 1; i <= ipVars_.maxIter; i++)
    {
      updatePinInfo();
      createSparseMatrix();

      // BiCGSTAB solver for initial place
      BiCGSTAB<SMatrix, IdentityPreconditioner> solver;
      solver.setMaxIterations(ipVars_.maxSolverIter);
      solver.compute(placeInstForceMatrixX_);
      instLocVecX_ = solver.solveWithGuess(fixedInstForceVecX_, instLocVecX_);
      errorX = solver.error();

      solver.compute(placeInstForceMatrixY_);
      instLocVecY_ = solver.solveWithGuess(fixedInstForceVecY_, instLocVecY_);
      errorY = solver.error();

      LOG_DEBUG("[InitialPlace] Iter {} CG Error: {} HPWL: {}", i, max(errorX, errorY), pb_->hpwl());
      updateCoordi();

      Plot::plot(pb_.get(), "./plot/cell", "ip_" + to_string(i));

      if (max(errorX, errorY) <= 1e-5 && i >= 5)
      {
        break;
      }
    }

    LOG_TRACE("Initial Placement End");
  }

  // starting point of initial place is center.
  void InitialPlace::placeInstsCenter()
  {
    const int centerX = pb_->dies().front()->coreCx();
    const int centerY = pb_->dies().front()->coreCy();

    for (auto &inst : pb_->insts())
    {
      if(!inst->isFixed())
        inst->setCenterLocation(centerX, centerY);
    }
  }

  void InitialPlace::setPlaceInstExtId()
  {
    // reset ExtId for all instances
    for (auto &inst : pb_->insts())
    {
      inst->setExtId(INT_MAX);
    }
    // set index only with place-able instances
    int id = 0;
    for (Instance* inst : pb_->insts())
    {
      if(!inst->isFixed())
      {
        inst->setExtId(id);
        id++;
      }
    }
  }

  void InitialPlace::updatePinInfo()
  {
    // reset all MinMax attributes
    for (auto &pin : pb_->pins())
    {
      pin->setMinPinX(false);
      pin->setMinPinY(false);
      pin->setMaxPinX(false);
      pin->setMaxPinY(false);
    }

    for (auto &net : pb_->nets())
    {
      Pin *pinMinX = nullptr, *pinMinY = nullptr;
      Pin *pinMaxX = nullptr, *pinMaxY = nullptr;
      int lx = INT_MAX, ly = INT_MAX;
      int ux = INT_MIN, uy = INT_MIN;

      // Mark B2B info on Pin structures
      for (auto &pin : net->pins())
      {
        if (lx > pin->cx())
        {
          if (pinMinX)
          {
            pinMinX->setMinPinX(false);
          }
          lx = pin->cx();
          pinMinX = pin;
          pinMinX->setMinPinX(true);
        }

        if (ux < pin->cx())
        {
          if (pinMaxX)
          {
            pinMaxX->setMaxPinX(false);
          }
          ux = pin->cx();
          pinMaxX = pin;
          pinMaxX->setMaxPinX(true);
        }

        if (ly > pin->cy())
        {
          if (pinMinY)
          {
            pinMinY->setMinPinY(false);
          }
          ly = pin->cy();
          pinMinY = pin;
          pinMinY->setMinPinY(true);
        }

        if (uy < pin->cy())
        {
          if (pinMaxY)
          {
            pinMaxY->setMaxPinY(false);
          }
          uy = pin->cy();
          pinMaxY = pin;
          pinMaxY->setMaxPinY(true);
        }
      }
    }
  }

  // solve placeInstForceMatrixX_ * xcg_x_ = xcg_b_ and placeInstForceMatrixY_ * ycg_x_ = ycg_b_ eq.
  void InitialPlace::createSparseMatrix()
  {
    instLocVecX_.resize(placeInstCnt_);
    fixedInstForceVecX_.resize(placeInstCnt_);
    instLocVecY_.resize(placeInstCnt_);
    fixedInstForceVecY_.resize(placeInstCnt_);

    placeInstForceMatrixX_.resize(placeInstCnt_, placeInstCnt_);
    placeInstForceMatrixY_.resize(placeInstCnt_, placeInstCnt_);

    //
    // listX and listY is a temporary vector that have tuples, (idx1, idx2, val)
    //
    // listX finally becomes placeInstForceMatrixX_
    // listY finally becomes placeInstForceMatrixY_
    //
    // The triplet vector is recommended usages
    // to fill in SparseMatrix from Eigen docs.
    //

    vector<T> listX, listY;
    listX.reserve(1000000);
    listY.reserve(1000000);

    // initialize vector
    for (Instance* inst : pb_->insts())
    {
      if(!inst->isFixed())
      {
        int idx = inst->extId();

        instLocVecX_(idx) = inst->cx();
        instLocVecY_(idx) = inst->cy();

        fixedInstForceVecX_(idx) = fixedInstForceVecY_(idx) = 0;
      }
    }

    // for each net
    for (auto &net : pb_->nets())
    {

      // skip for small nets.
      if (net->pins().size() <= 1)
      {
        continue;
      }

      // escape long time cals on huge fanout.
      //
      if (net->pins().size() >= ipVars_.maxFanout)
      {
        continue;
      }

      float netWeight = ipVars_.netWeightScale / (net->pins().size() - 1);
      // cout << "net: " << net.net()->getConstName() << endl;

      // foreach two pins in single nets.
      for (auto &pin1 : net->pins())
      {
        int pinIdx1 = &pin1 - &(net->pins()[0]);
        for (auto &pin2 : net->pins())
        {
          int pinIdx2 = &pin2 - &(net->pins()[0]);

          //
          // will compare two pins "only once."
          //
          if (pinIdx1 < pinIdx2)
          {
            break;
          }

          // no need to fill in when instance is same
          if (pin1->instance() == pin2->instance())
          {
            continue;
          }

          // B2B modeling on min/maxX pins.
          if (pin1->isMinPinX() || pin1->isMaxPinX() ||
              pin2->isMinPinX() || pin2->isMaxPinX())
          {
            int diffX = abs(pin1->cx() - pin2->cx());
            float weightX = 0;
            if (diffX > ipVars_.minDiffLength)
            {
              weightX = netWeight / diffX;
            }
            else
            {
              weightX = netWeight / ipVars_.minDiffLength;
            }
            // cout << weightX << endl;

            // both pin cames from instance
            if (!pin1->instance()->isFixed() && !pin2->instance()->isFixed())
            {
              const int inst1 = pin1->instance()->extId();
              const int inst2 = pin2->instance()->extId();
              // cout << "inst: " << inst1 << " " << inst2 << endl;

              listX.push_back(T(inst1, inst1, weightX));
              listX.push_back(T(inst2, inst2, weightX));

              listX.push_back(T(inst1, inst2, -weightX));
              listX.push_back(T(inst2, inst1, -weightX));

              // cout << pin1->cx() << " "
              //   << pin1->instance()->cx() << endl;
              fixedInstForceVecX_(inst1) +=
                  -weightX * ((pin1->cx() - pin1->instance()->cx()) -
                              (pin2->cx() - pin2->instance()->cx()));

              fixedInstForceVecX_(inst2) +=
                  -weightX * ((pin2->cx() - pin2->instance()->cx()) -
                              (pin1->cx() - pin1->instance()->cx()));
            }
            // pin1 from IO port / pin2 from Instance
            else if (pin1->instance()->isFixed() && !pin2->instance()->isFixed())
            {
              const int inst2 = pin2->instance()->extId();
              // cout << "inst2: " << inst2 << endl;
              listX.push_back(T(inst2, inst2, weightX));
              fixedInstForceVecX_(inst2) += weightX *
                                            (pin1->cx() -
                                             (pin2->cx() - pin2->instance()->cx()));
            }
            // pin1 from Instance / pin2 from IO port
            else if (!pin1->instance()->isFixed() && pin2->instance()->isFixed())
            {
              const int inst1 = pin1->instance()->extId();
              // cout << "inst1: " << inst1 << endl;
              listX.push_back(T(inst1, inst1, weightX));
              fixedInstForceVecX_(inst1) += weightX *
                                            (pin2->cx() -
                                             (pin1->cx() - pin1->instance()->cx()));
            }
          }

          // B2B modeling on min/maxY pins.
          if (pin1->isMinPinY() || pin1->isMaxPinY() ||
              pin2->isMinPinY() || pin2->isMaxPinY())
          {

            int diffY = abs(pin1->cy() - pin2->cy());
            float weightY = 0;
            if (diffY > ipVars_.minDiffLength)
            {
              weightY = netWeight / diffY;
            }
            else
            {
              weightY = netWeight / ipVars_.minDiffLength;
            }

            // both pin cames from instance
            if (!pin1->instance()->isFixed() && !pin2->instance()->isFixed())
            {
              const int inst1 = pin1->instance()->extId();
              const int inst2 = pin2->instance()->extId();

              listY.push_back(T(inst1, inst1, weightY));
              listY.push_back(T(inst2, inst2, weightY));

              listY.push_back(T(inst1, inst2, -weightY));
              listY.push_back(T(inst2, inst1, -weightY));

              fixedInstForceVecY_(inst1) +=
                  -weightY * ((pin1->cy() - pin1->instance()->cy()) -
                              (pin2->cy() - pin2->instance()->cy()));

              fixedInstForceVecY_(inst2) +=
                  -weightY * ((pin2->cy() - pin2->instance()->cy()) -
                              (pin1->cy() - pin1->instance()->cy()));
            }
            // pin1 from IO port / pin2 from Instance
            else if (pin1->instance()->isFixed() && !pin2->instance()->isFixed())
            {
              const int inst2 = pin2->instance()->extId();
              listY.push_back(T(inst2, inst2, weightY));
              fixedInstForceVecY_(inst2) += weightY *
                                            (pin1->cy() -
                                             (pin2->cy() - pin2->instance()->cy()));
            }
            // pin1 from Instance / pin2 from IO port
            else if (!pin1->instance()->isFixed() && pin2->instance()->isFixed())
            {
              const int inst1 = pin1->instance()->extId();
              listY.push_back(T(inst1, inst1, weightY));
              fixedInstForceVecY_(inst1) += weightY *
                                            (pin2->cy() -
                                             (pin1->cy() - pin1->instance()->cy()));
            }
          }
        }
      }
    }

    placeInstForceMatrixX_.setFromTriplets(listX.begin(), listX.end());
    placeInstForceMatrixY_.setFromTriplets(listY.begin(), listY.end());
  }

  void InitialPlace::updateCoordi()
  {
    for (Instance* inst : pb_->insts())
    {
      if(!inst->isFixed())
      {
        int idx = inst->extId();
        inst->setCenterLocation(instLocVecX_(idx), instLocVecY_(idx));
      }
    }
  }

}
