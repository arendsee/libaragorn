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

class Gene {
  public:
    int ps;
    int nbase;
    int start;
    int stop;
    int astem1;
    int astem2;
    int spacer1;
    int spacer2;
    int dstem;
    int dloop;
    int cstem;
    int cloop;
    int intron;
    int nintron;
    int anticodon;
    int var;
    int tstem;
    int tloop;
    double energy;
    int tps;
    int tpe;
    int asst;

  Gene () {
    ps        = 0;
    nbase     = 0;
    start     = 0L;
    stop      = 0L;
    astem1    = 7;
    astem2    = 7;
    spacer1   = 2;
    spacer2   = 1;
    dstem     = 3;
    dloop     = 9;
    cstem     = 5;
    cloop     = 7;
    intron    = 0;
    nintron   = 0;
    anticodon = 0;
    var       = 15;
    tstem     = 5;
    tloop     = 7;
    energy    = 0.0;
    tps       = -1; // just for tmRNA
    tpe       = -1; // just for tmRNA
    asst      = -1;
  }
};

class TrnaLoop {
  public:
    int pos;
    int stem; // length of the T-stem (4 or 5)
    int loop;
    double energy;
};

class TrnaAstem {
  public:
    int pos;
    double energy;
};

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

// find all tstems in the sequence
std::vector<TrnaLoop> find_tstems(const std::vector<int>& s, int loffset, int roffset, double ttscanthresh, double ttarmthresh);

std::vector<TrnaAstem> find_astem5(const std::vector<int>& seq, int si, int sl, int astem3, int n3, double tascanthresh, double tastemthresh);

tRNA make_trna(Gene &g);

TrnaLoop make_trna_loop(int pos, int loop, int stem, double energy);

TrnaAstem make_astem(int pos, double energy);

#endif
