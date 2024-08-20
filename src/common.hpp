#ifndef __ARAGORN_COMMON_HPP__
#define __ARAGORN_COMMON_HPP__

#include <vector>
#include <string>
#include "aragorn.hpp"

#define Adenine         0
#define Cytosine        1
#define Guanine         2
#define Thymine         3
#define AMBIG           4
#define NOBASE          5

#define ND          100
#define ATBOND      2.5
#define GCBOND      3.0

#define ASTEM2_EXT      9

/* Basepair matching matrices */

static int bp[6][6] =
 { { 0,0,0,1,1,0 },
   { 0,0,1,0,1,0 },
   { 0,1,0,1,1,0 },
   { 1,0,1,0,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

static int vbp[6][6] =
 { { 0,0,1,4,4,0 },
   { 0,0,4,0,4,0 },
   { 1,4,0,2,4,0 },
   { 4,0,2,0,4,0 },
   { 4,4,4,4,4,0 },
   { 0,0,0,0,0,0 } };

std::vector<int> dna2int (const std::string& dna);

std::vector<tRNA> best_hit(std::vector<tRNA> hits);

#endif
