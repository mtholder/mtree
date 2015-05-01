#if !defined(__PRUNE_STRUCT_H__)
#define __PRUNE_STRUCT_H__
#include <cassert>
#include <stack>
#include <climits>
#include <iostream>
#include <vector>
#include "mt_log.h"
#include "mt_part_data.h"
#include "mt_part_model_info.h"
namespace mt {

// taken from local vars of pllNewviewIterative in newviewGenericSpecial.c
class PruneStruct {
    double *x1_start;
    double *x2_start;
    double *x3_start;
    double *left;
    double *right;
    #if (defined(__SSE3) || defined(__AVX))
        double *x1_gapColumn;
        double *x2_gapColumn;
        double *x3_gapColumn;
        unsigned int * x1_gap;
        unsigned int * x2_gap;
        unsigned int * x3_gap;
        std::size_t gapOffset;
    #endif
    unsigned char * tipX1;
    unsigned char * tipX2;
    const double*rateCategories;
    double*x1_ascColumn;
    double*x2_ascColumn;
    double*x3_ascColumn;
    int categories;
    int scalerIncrement;
    std::size_t rateHet;
    std::size_t ascWidth;
    std::size_t states;
    std::size_t availableLength;
    std::size_t requiredLength;
    const int * wgt; // integer weight vector with pattern compression weights

    PruneStruct(const PartData & pd,
                const SlotIndices & slots)
        :x1_start(nullptr),
        x2_start(nullptr),
        x3_start(pd.GetCLArray(p_slot)),
        left(nullptr),
        right(nullptr),
#       if (defined(__SSE3) || defined(__AVX))
            x1_gapColumn(nullptr),
            x2_gapColumn(nullptr),
            x3_gapColumn(nullptr),
            x1_gap(nullptr),
            x2_gap(nullptr),
            x3_gap(nullptr),
            gapOffset(0),
#       endif
        tipX1(nullptr),
        tipX2(nullptr),
        rateCategories(psi.GetRateCategories()),
        x1_ascColumn(nullptr),
        x2_ascColumn(nullptr),
        x3_ascColumn(nullptr),
        categories(psi.GetNumRateHetCategories()),
        scalerIncrement(0),
        rateHet(psi.GetNumRateHetCategories()),
        ascWidth(psi.GetAscWidth()),
        states(pd->GetNumStates()),
            // get the length of the current likelihood array stored at node p. This is 
            // important mainly for the SEV-based memory saving option described in here:
            // F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees".
            // So pr->partitionData[model]->xSpaceVector[i] provides the length of the allocated conditional array of partition model
            // and node i 
        availableLength(pd.GetXSpaceVector(slots.p)),
        requiredLength(0),
        wgt(pd.GetIntWeigths()) {
            // if we are not trying to save memory the space required to store an inner likelihood array 
            // is the number of sites in the partition times the number of states of the data type in the partition 
            // times the number of discrete GAMMA rates (1 for CAT essentially) times 8 bytes
            requiredLength  =  virtual_width( width ) * rateHet * states * sizeof(double);
            //                   printf( "req: %d %d %d %d\n", requiredLength, width, virtual_width(width), model );
            // Initially, even when not using memory saving no space is allocated for inner likelihood arrats hence 
            //  availableLength will be zero at the very first time we traverse the tree.
            //  Hence we need to allocate something here 
            if(requiredLength != availableLength) {
                pd.allocXVector(ps);
            }
        }

};


    void allocXVector(PruneStruct & ps) {
        /* if there is a vector of incorrect length assigned here i.e., x3 != NULL we must free  it first */
        if(x3_start)
            rax_free(x3_start);
        /* allocate memory: note that here we use a byte-boundary aligned malloc, because we need the vectors
        to be aligned at 16 BYTE (SSE3) or 32 BYTE (AVX) boundaries! */
        rax_posix_memalign ((void **)&ps.x3_start, PLL_BYTE_ALIGNMENT, ps.requiredLength);              
        /* update the data structures for consistent bookkeeping */
        internalCLVector[p_slot]      = ps.x3_start;
        xSpaceVector[p_slot] = ps.requiredLength;
    }

} //namespace
#endif
