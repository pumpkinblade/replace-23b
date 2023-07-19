#include "randomPartitioner.h"
#include "placerBase.h"
#include "placer23b.h"

namespace replace {

PlacerBase RandomPartitioner::doPartition(){
    PlacerBase pb;

    // make some space for inst, net and pin.
    pb.instStor_.reserve(p23_.instStor_.size());
    pb.netStor_.reserve(p23_.netStor_.size());
    pb.pinStor_.reverse(p23_.pinStor_.size());

    // create top, botton, terminal dies with 
    // same shape and row params of placer23
    for(int i=0; i<3; i++){
        pb.dieStor_.emplace_back();
        Die& dieref = pb.dieStor_.back();
        switch (i)
        {
        case 0: dieref.copyDieBoxAndRowParamsFrom(p23_.topDie());
            break;
        case 1: dieref.copyDieBoxAndRowParamsFrom(p23_.boottomDie());
            break;
        // for terminal die, only set die box.
        case 2: dieref.setDieBox(p23_.topDie()); break;
        default: break;
        }
        pb.dies_.push_back(&dieref);
    }

    // randomly assign instance to top or bottom. Set instance shape by 
    // technology of the die it assigned to.
    srand(10086);
    for(auto inst : p23_.instStor_){
        float roll = (float) rand()/RAND_MAX;
        bool assignToTop = roll >= 0.5 ? true : false;
        if(assignToTop){ pb.dies()[0]->addInstance(&inst); }
        else{ pb.dies()[1]->addInstance(&inst); }
    }

    // TODO: create terminal on terminal die if a net is seperated to two dies

    // derive pins_, nets_, insts_ of pb
    pb.deriveIPNs()

    return pb;
}
}
