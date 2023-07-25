#ifndef __REPLACE_MACRO_LEGALIZER__
#define __REPLACE_MACRO_LEGALIZER__

#include <memory>
#include <vector>
#include <cassert>

namespace replace{
    class Instance;
    class Die;
    class PlacerBase;
    class NesterovBase;
    class GCell;

    class MacroLegalizerVars
    {
        public:
            float t, r;


            MacroLegalizerVars();
    };

    class macroLegalizer
    {
        public:
            macroLegalizer(/* args */);
            ~macroLegalizer();
            void doMacroLegalization();
            int overlapArea(GCell* , GCell*);
            int getCellMacroOverlap();
            int getMacrosOverlap();

        private:
        /* data */
            std::shared_ptr<PlacerBase> pb_;
            std::shared_ptr<NesterovBase> nb_;
            MacroLegalizerVars vars_;
    };
    
    

}

#endif