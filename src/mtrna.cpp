#include "aragorn.hpp"
#include "common.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>

#define MAMMAL_MT  2

#define TERM            -1
#define INACTIVE        2.0e+35

#define ASTEM2_EXT      9
#define MAXINTRONLEN    3000  // FIXME: introns should not be present in mitochondrial tRNAs
#define MINCTRNALEN     62
#define MAXCTRNALEN     110
#define MINTRNALEN      (MINCTRNALEN + 1)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM2_EXT)
#define MAXETRNALEN     (MAXTRNALEN + MAXINTRONLEN)

#define mtNA        1500
#define mtND        150 
#define mtNTH       3000 
#define mtNTM       3
#define mtNCDS      200
#define mtNCDSCODON 6000
#define mtGCBOND    0.0
#define mtATBOND    -0.5
#define mtGTBOND    -1.2
#define mtTTBOND    -2.9
#define mtGGBOND    -3.0
#define mtGABOND    -3.0
#define mtNOBOND    -3.0
#define mtBONDSTAB  1.5
#define mtABONDSTAB 2.0
#define mtTSTTSTAB   -2.5
#define mtTERMSTAB  0.01
#define mtSENDSTAB  0.01
#define mtNSTAB     0.1
#define mt3MMSTAB   1.0
#define mtGCPENALTY 0.8
#define mtGCPENALTYD 2.0
#define mt_DRLmaxlength 16
#define mt_TVRLmaxlength 18
#define mtNCLM 3

typedef struct { int *pos;
                 int stem;
                 int loop;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_loop;

typedef struct { int *pos;
                 int *looppos;
                 int *end;
                 int stem;
                 int loop;
                 int arm;
                 int anticodon;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_cloop;

typedef struct { int *pos;
                 int stem;
                 int loop;
		         int *end;
                 unsigned int bondtype;
                 double energy;
                 double stem_energy; } mt_trna_tloop;

typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 unsigned int bondtype;
                 double energy; } mt_trna_astem;

typedef struct { int *pos1;
                 int *pos2;
                 int comp; } mt_cds;

typedef struct { char name[100];
                 int seq[MAXTRNALEN+1];
                 int eseq[MAXETRNALEN+1];
                 int *ps;
                 int nbase;
                 int comp;
                 long start;
                 long stop;
                 int astem1;
                 int astem2;
                 int aatail;
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
                 int varbp;
                 int tstem;
                 int tloop;
                 int genetype;
                 double energy;
                 int asst;
                 int tps;
                 int tpe;
                 int annotation;
                 int annosc;   } gene;

int lbp[3][6][6] =
 { { { 0,0,1,1,1,0 },
     { 0,0,1,0,1,0 },
     { 1,1,0,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } },
   { { 0,0,0,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,0,1,1,0 },
     { 1,0,1,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } },
   { { 0,0,0,1,1,0 },
     { 0,0,1,0,1,0 },
     { 0,1,0,0,1,0 },
     { 1,0,0,0,1,0 },
     { 1,1,1,1,1,0 },
     { 0,0,0,0,0,0 } } };

int wbp[6][6] =
 { { 0,0,0,2,2,0 },
   { 0,0,2,0,2,0 },
   { 0,2,0,1,2,0 },
   { 2,0,1,0,2,0 },
   { 2,2,2,2,2,0 },
   { 0,0,0,0,0,0 } };

int wcbp[6][6] = { { 0,0,0,1,1,0 },
                   { 0,0,1,0,1,0 },
                   { 0,1,0,0,1,0 },
                   { 1,0,0,0,1,0 },
                   { 1,1,1,1,1,0 },
                   { 0,0,0,0,0,0 } };

int gc[6][6] = { { 0,0,0,0,0,0 },
                 { 0,0,1,0,1,0 },
                 { 0,1,0,0,1,0 },
                 { 0,0,0,0,0,0 },
                 { 0,1,1,0,1,0 },
                 { 0,0,0,0,0,0 } };

int gt[6][6] = { { 0,0,0,0,0,0 },
                 { 0,0,0,0,0,0 },
                 { 0,0,0,1,1,0 },
                 { 0,0,1,0,1,0 },
                 { 0,0,1,1,1,0 },
                 { 0,0,0,0,0,0 } };

int at[6][6] = { { 0,0,0,1,1,0 },
                 { 0,0,0,0,0,0 },
                 { 0,0,0,0,0,0 },
                 { 1,0,0,0,0,0 },
                 { 1,0,0,0,1,0 },
                 { 0,0,0,0,0,0 } };

int tt[6][6] = { { 0,0,0,0,0,0 },
                 { 0,0,0,0,0,0 },
                 { 0,0,0,0,0,0 },
                 { 0,0,0,1,1,0 },
                 { 0,0,0,1,1,0 },
                 { 0,0,0,0,0,0 } };

int stemterm[6][6] = { { 0,0,1,0,1,0 },
                       { 0,0,0,0,0,0 },
                       { 1,0,0,0,1,0 },
                       { 0,0,0,1,1,0 },
                       { 1,0,1,1,1,0 },
                       { 0,0,0,0,0,0 } };

int aastemterm[6][6] =
 { { 1,0,1,0,1,0 },
   { 0,0,0,0,0,0 },
   { 1,0,0,0,1,0 },
   { 0,0,0,1,1,0 },
   { 1,0,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int ggstemterm[6][6] =
 { { 0,0,1,0,1,0 },
   { 0,0,0,0,0,0 },
   { 1,0,1,0,1,0 },
   { 0,0,0,1,1,0 },
   { 1,0,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int assymst[6][6] = { { 0,0,0,0,0,0 },
                      { 0,0,0,0,0,0 },
                      { 1,0,0,0,1,0 },
                      { 0,0,0,1,1,0 },
                      { 1,0,0,1,1,0 },
                      { 0,0,0,0,0,0 } };

int assymat[6][6] = { { 0,0,0,1,1,0 },
                      { 0,0,0,0,0,0 },
                      { 0,0,0,0,0,0 },
                      { 0,0,0,0,0,0 },
                      { 0,0,0,1,1,0 },
                      { 0,0,0,0,0,0 } };


int stackbp[6][6] = { { 0,0,0,1,1,0 },
                      { 0,0,1,0,1,0 },
                      { 0,1,0,1,1,0 },
                      { 1,0,1,1,1,0 },
                      { 1,1,1,1,1,0 },
                      { 0,0,0,0,0,0 } };

int ggstackbp[6][6] =
 { { 0,0,0,1,1,0 },
   { 0,0,1,0,1,0 },
   { 0,1,1,1,1,0 },
   { 1,0,1,1,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };


int ggbp[6][6] =
 { { 0,0,0,1,1,0 },
   { 0,0,1,0,1,0 },
   { 0,1,1,1,1,0 },
   { 1,0,1,0,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int gabp[6][6] =
 { { 0,0,1,1,1,0 },
   { 0,0,1,0,1,0 },
   { 1,1,0,1,1,0 },
   { 1,0,1,0,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int assymagbp[6][6] =
 { { 0,0,1,1,1,0 },
   { 0,0,1,0,1,0 },
   { 0,1,0,1,1,0 },
   { 1,0,1,0,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int stembp[6][6] =
 { { 0,0,1,1,1,0 },
   { 0,0,1,0,1,0 },
   { 1,1,0,1,1,0 },
   { 1,0,1,1,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int ggstembp[6][6] =
 { { 0,0,1,1,1,0 },
   { 0,0,1,0,1,0 },
   { 1,1,1,1,1,0 },
   { 1,0,1,1,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int gastembp[6][6] =
 { { 1,0,1,1,1,0 },
   { 0,0,1,0,1,0 },
   { 1,1,1,1,1,0 },
   { 1,0,1,1,1,0 },
   { 1,1,1,1,1,0 },
   { 0,0,0,0,0,0 } };

int tandemid[mtNTM][4] =
 { { 3,2,2,3 },
   { 2,3,3,2 },
   { 3,3,3,3 } };

double tandem_em[mtNTM] = { -0.5,-0.5,2.0 };

double send_em[6][6] =
 { { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,0.5*mtSENDSTAB,0.0,0.5*mtSENDSTAB,0.0 },
   { 0.0,0.5*mtSENDSTAB,0.0,mtSENDSTAB,mtSENDSTAB,0.0 },
   { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
   { 0.0,0.5*mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 } };

double ssend_em[6][6] =
 { { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
   { 0.0,mtSENDSTAB,0.0,mtSENDSTAB,mtSENDSTAB,0.0 },
   { 0.0,0.0,mtSENDSTAB,0.0,mtSENDSTAB,0.0 },
   { 0.0,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,mtSENDSTAB,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 } };

int neighbour_map[6][6] =
 { { 0,0,1,0,1,0 },
   { 0,0,0,0,0,0 },
   { 1,0,0,0,1,0 },
   { 0,0,0,1,1,0 },
   { 1,0,1,1,1,0 },
   { 0,0,0,0,0,0 } };

double neighbour_em[2][6][6] = {
 { { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 } },

 { { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,0.0,mtNSTAB,0.0,mtNSTAB,0.0 },
   { 0.0,mtNSTAB,0.0,0.0,mtNSTAB,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 },
   { 0.0,mtNSTAB,mtNSTAB,0.0,mtNSTAB,0.0 },
   { 0.0,0.0,0.0,0.0,0.0,0.0 } } };

unsigned int btmap[6][6] =
 { { 0x10000,0x10000,0x1000,0x10,0x00000,0x10000 },
   { 0x10000,0x10000,0x1,0x10000,0x00000,0x10000 },
   { 0x1000,0x1,0x10000,0x100,0x00000,0x10000 },
   { 0x10,0x10000,0x100,0x1000,0x00000,0x10000 },
   { 0x00000,0x00000,0x00000,0x00000,0x00000,0x10000 },
   { 0x10000,0x10000,0x10000,0x10000,0x10000,0x10000 } };

int mt_discrim[3][64][6] =
 /* metazoan mt */
 {{{ 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 0,0,0,0,0,0 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 0,0,0,0,0,0 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 }},
/* standard */
  {{ 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },

   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 },  { 1,1,1,1,1,1 }},
 /* mammal mt */
  {{ 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },
   { 0,0,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
   { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,1,0,0,1,1 },  { 0,0,1,0,1,1 },

   { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,0,1,1 },
   { 0,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,0,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
   { 0,0,0,0,0,0 },  { 1,0,1,1,1,1 },  { 0,0,1,0,1,1 },  { 1,1,1,1,1,1 },

   { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,0,1,1 },
   { 0,1,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },
   { 1,1,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,1,0,0,1,1 },  { 0,0,1,0,1,1 },

   { 1,1,0,1,1,1 },  { 1,0,0,0,1,1 },  { 1,0,1,1,1,1 },  { 1,0,0,0,1,1 },
   { 1,0,1,0,1,1 },  { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 0,0,0,0,0,0 },
   { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,0,1,0,1,1 },  { 1,1,1,1,1,1 },
   { 1,0,0,0,1,1 },  { 1,1,1,1,1,1 },  { 1,0,1,0,1,1 },  { 1,1,1,1,1,1 }}};

double par_mttarmthresh = -7.9;
double par_mtdtthresh = 91.5;
double par_mttthresh = 83.5;
double par_mtdthresh = 85.0;
int par_discrim = 0; // metazoan table
int par_cloop7 = 0;
int par_mtxdetect = 1; // see -x option which sets this to 0
int par_loffset = 0;
int par_roffset = 0;
int par_tvloop = 1;
// int par_aataildiv = 0; // 1 allows some divergence from the 3' amino-acyl acceptor sequence NCCA (-jr4 option)


double vloop_stability(int *sb, int var, int *varbp)
{ int e,stem,vstem,loop,*sn,*sen,*pos1,*pos2,*se,*sc,*sd,*sf,*s;
  unsigned int c,cn,m;
  static unsigned int A[6] = { 0,0,0x100,0x400,0,0 };
  static unsigned int C[6] = { 0,0,0x400,0,0,0 };
  static unsigned int G[6] = { 0x100,0x400,0,0x200,0,0 };
  static unsigned int T[6] = { 0x400,0,0x200,0,0,0 };
  static unsigned int te[6] = { 0,0,0,0,0,0 };
  e = 0;
  sc = sb + 3;   
  se = sb + var - 2; 
  sf = se - 2;
  te[0] = A[*se];
  te[1] = C[*se];
  te[2] = G[*se];
  te[3] = T[*se];
  while (--se > sf)
   { te[0] = (te[0] >> 4) | A[*se];
     te[1] = (te[1] >> 4) | C[*se];
     te[2] = (te[2] >> 4) | G[*se];
     te[3] = (te[3] >> 4) | T[*se]; }
  while (se >= sc)
   { te[0] = ((te[0] >> 4) | A[*se]);
     te[1] = ((te[1] >> 4) | C[*se]);
     te[2] = ((te[2] >> 4) | G[*se]);
     te[3] = ((te[3] >> 4) | T[*se]);
     s = se - 5;
     sd = se - 7;
     m = te[*s];
     while (--s > sd) m = (m >> 4) + te[*s];
     while (s >= sb)
       {  m = (m >> 4) + te[*s];
          c = m & 0xf;
          if (c >= 9)
           { stem = 3;
             loop = (int)(se - s) - 3;
             sen = se;
             sn = s + 2;
             while (loop >= 6)
              { if ((cn = vbp[sen[-1]][sn[1]]) <= 0) break;
                c += cn;
                stem++;
                loop -= 2;
                sen--;
                sn++; }
             if (c > e)
              { e = c;
                pos1 = s;
                pos2 = sen;
                vstem = stem; }}
          s--; }
      se--; }
  if (e > 0)
   { *varbp = (((int)(pos1-sb))<<10) + (((int)(pos2-sb))<<5) + vstem;
     return((double)(3*(vstem - 4))); }
  else
   { *varbp = 0;
     return(-12.0); }}

tRNA make_trna(gene &g){
    tRNA h;
    h.start = g.start;
    h.stop = g.stop;
    h.astem1 = g.astem1;
    h.astem2 = g.astem2;
    h.spacer1 = g.spacer1;
    h.spacer2 = g.spacer2;
    h.dstem = g.dstem;
    h.dloop = g.dloop;
    h.cstem = g.cstem;
    h.cloop = g.cloop;
    h.intron_start = g.intron;
    h.intron_length = g.nintron;
    h.anticodon = g.anticodon;
    h.var = g.var;
    h.tstem = g.tstem;
    h.tloop = g.tloop;
    h.score = g.energy;
    return h;
}

std::vector<tRNA> predict_mtrnas(std::string &dna) {

  std::vector<int> seq_vec = dna2int(dna);
  int* seq = &seq_vec[0];
  int lseq = seq_vec.size();

  std::vector<tRNA> gs;

  int nah,ndh,nch,nth,ncdsh,h,i,j,k,n,p,y,av,gcc,cgcc,catc,athresh;
  int igc,nbase,b8,b9,b48,b57,nc,na,nt,nti,nd,ndi,dposmap[32];
  int dl,tl,astem8,ti,di;
  int astem,asteme,cloop,dloop,tloop,tc;
  int carm,cstem,darm,dstem,tarm,tstem,var,varbp,spacer1,spacer2,anticodon;
  int dstemmotif,cloop7,mtxdetect,incds;
  int *s,*sl,*s1,*s2,*s4,*sa,*sb,*sc,*se,*sf,*sg,*si;
  int *slm,*slm1,*sle,*slb,*sld,*sge;
  int *dpos,*cpos,*cend,*tpos,*tend,*apos1,*apos2,*aend1,*aend2;
  int *clooppos,*cloopend;
  unsigned int bondtype,abondtype,mabondtype,acbondtype,cbondtype;
  unsigned int agcat,cgcat,tgcat,dbondtype,dtbondtype,tbondtype;
  unsigned int r,ct[6],cm,cv,q,tendmap[63];
  double gcv,e,ec,ea,eas,ed,et,ev,energy,stem_energy;
  double tarmthresh,tthresh,dthresh,dtthresh,thresh;
  mt_trna_cloop chit[6];
  static mt_trna_loop dhit[mtND+1];
  static mt_trna_tloop thit[mtNTH+1];
  static mt_trna_astem ahit[mtNA+1];
  static mt_cds cdshit[mtNCDS];
  static gene te =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,2,1,4,7,5,7,0,0,0,5,0,5,7,
     0,0.0,0,0,0 };
  static int cAI[6] = { 8,0,0,0,8,0 };
  static int cRI[6] = { 8,0,4,0,8,0 };
  static int cTI[6] = { 0,0,0,16,16,0 };
  static int cYI[6] = { 0,8,0,4,8,0 };
  static int AI[6] = { 1,0,0,0,1,0 };
  static int CI[6] = { 0,1,0,0,1,0 };
  static int GI[6] = { 0,0,1,0,1,0 };
  static int TI[6] = { 0,0,0,1,1,0 };
  static int RI[6] = { 1,0,1,0,1,0 };
  static int YI[6] = { 0,1,0,1,1,0 };
  static unsigned int tem[6] = { 0,0,0,0,0,0 };
  static unsigned int At[6] = { 0,0,0,1,1,0 };
  static unsigned int Ct[6] = { 0,0,1,0,1,0 };
  static unsigned int Gt[6] = { 0,1,0,1,1,0 };
  static unsigned int Tt[6] = { 1,0,1,0,1,0 };
  static unsigned int cAt[6] = { 0,0,0,2,2,0 };
  static unsigned int cCt[6] = { 0,0,2,0,2,0 };
  static unsigned int cGt[6] = { 0,2,0,1,2,0 };
  static unsigned int cTt[6] = { 2,0,1,0,2,0 };
  static unsigned int aAt[6] = { 0,0,1,2,2,0 };
  static unsigned int aCt[6] = { 0,0,2,0,2,0 };
  static unsigned int aGt[6] = { 1,2,0,1,2,0 };
  static unsigned int aTt[6] = { 2,0,1,1,2,0 };
  static unsigned int dAt[6] = { 0,0,1,2,2,0 };
  static unsigned int dCt[6] = { 0,0,2,0,2,0 };
  static unsigned int dGt[6] = { 1,2,0,2,2,0 };
  static unsigned int dTt[6] = { 2,0,2,1,2,0 };
  static unsigned int clmotif[mtNCLM] = { 0x1321300,0x3321300,0x1323002 };
  static int dloopi[mt_DRLmaxlength+1][4] =
   { { -1 }, { -1 }, { -1 }, { -1 }, { -1 }, { -1 }, { -1 },
     { 0,2,-1 }, { 0,2,-1 }, { 0,2,3,-1 }, { 0,3,-1 }, { 0,3,-1 },
     { 0,3,4,-1 }, { 0,4,-1 }, { 0,5,-1 }, { 0,5,6,-1 }, { 0,5,6,-1 } };
  static int tloopa[12][4] =
   { { -1 }, { -1 }, { -1 }, { 0,1,-1 }, { 0,2,1,-1 }, { 4,3,2,-1 },
     { 4,3,-1 }, { 4,3,-1 }, { 4,3,-1 }, { 5,4,3,-1 }, { 5,4,-1 }, { 5,-1 } };
  static double dA[6] = { 1.0,0.0,0.0,0.0,1.0,0.0 };
  static double dT[6] = { 0.0,0.0,0.0,1.0,1.0,0.0 };
  static double C[6] = { 0.0,1.0,0.0,0.0,1.0,0.0 };
  static double G[6] = { 0.0,0.0,1.0,0.0,1.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,1.0,1.0,0.0 };
  static double AX[6] = { 0.0,-1.0,-1.0,-1.0,0.0,-1.0 };
  static double AX37[6] = { 0.0,-4.0,-1.0,-4.0,0.0,-4.0 };
  static double AXX[6] = { 0.0,-3.0,-1.5,-3.0,0.0,-3.0 };
  static double AXX37[6] = { 0.0,-4.0,-4.0,-4.0,0.0,-4.0 };
  static double AX7[6] = { 0.0,-0.7,-0.7,-0.7,0.0,-0.7 };
  static double CX[6] = { -2.0,0.0,-2.0,-1.0,0.0,-2.0 };
  static double CX7[6] = { -0.7,0.0,-0.7,-0.7,0.0,-0.7 };
  static double TX[6] = { -1.0,-1.0,-1.0,0.0,0.0,-1.0 };
  static double tC[6] = { 0.0,0.01,0.0,0.0,0.01,0.0 };
  static double tG[6] = { 0.0,0.0,0.01,0.0,0.01,0.0 };
  static double tT[6] = { 0.0,0.0,0.0,0.01,0.01,0.0 };
  static double cA[6] = { 0.8,0.0,0.0,0.0,0.8,0.0 };
  static double cR[6] = { 0.8,-2.0,0.8,-0.8,0.8,-0.8 };
  static double cT[6] = { -0.8,0.0,-0.8,2.6,2.6,-0.8 };
  static double cY[6] = { -0.8,0.8,-0.8,0.8,0.8,-0.8 };
  static double loop_stab[41] =
  { 10.0,2.0,1.0,0.4,0.3,0.2,0.1,0.0,0.1,0.2,0.3,0.4,0.5,1.6,1.7,1.8,
    1.9,2.0,2.1,2.2,2.3,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,
    5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7 };
  static double bem[6][6] =
   { {  mtNOBOND, mtNOBOND, mtGABOND, mtATBOND, mtATBOND, mtNOBOND },
     {  mtNOBOND, mtNOBOND, mtGCBOND, mtNOBOND, mtGCBOND, mtNOBOND },
     {  mtGABOND, mtGCBOND, mtGGBOND, mtGTBOND, mtGCBOND, mtNOBOND },
     {  mtATBOND, mtNOBOND, mtGTBOND, mtTTBOND, mtATBOND, mtNOBOND },
     {  mtATBOND, mtGCBOND, mtGCBOND, mtATBOND, mtGCBOND, mtNOBOND },
     {  mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND, mtNOBOND } };
  static double hbem[5][5] =
   { {  0.0,0.0,0.0,mtBONDSTAB+0.5*mtATBOND,mtBONDSTAB+0.5*mtATBOND },
     {  0.0,0.0,mtBONDSTAB+0.5*mtGCBOND,0.0,mtBONDSTAB+0.5*mtGCBOND },
     {  0.0,mtBONDSTAB+0.5*mtGCBOND,0.0,mtBONDSTAB+0.5*mtGTBOND,
                                            mtBONDSTAB+0.5*mtGCBOND },
     {  mtBONDSTAB+0.5*mtATBOND,0.0,mtBONDSTAB+0.5*mtGTBOND,0.0,
                                            mtBONDSTAB+0.5*mtATBOND },
     {  mtBONDSTAB+0.5*mtATBOND,mtBONDSTAB+0.5*mtGCBOND,
                                            mtBONDSTAB+0.5*mtGCBOND,
                                            mtBONDSTAB+0.5*mtATBOND,
                                            mtBONDSTAB+0.5*mtGCBOND } };
  tarmthresh = par_mttarmthresh;
  tthresh = par_mttthresh;
  dthresh = par_mtdthresh;
  dtthresh = par_mtdtthresh;
  cloop7 = par_cloop7;
  mtxdetect = par_mtxdetect;
  ncdsh = 0;

  sc = seq + par_loffset;
  sl = seq + lseq - par_roffset;
  h = sc[16];
  p = sc[15];
  j = sc[14];
  k = sc[13];
  n = sc[12];
  y = sc[11];
  ct[0] = cAt[h] | (cAt[p]<<4) | (cAt[j]<<8) | (cAt[k]<<12) | (cAt[n]<<16) | (cAt[y]<<20);
  ct[1] = cCt[h] | (cCt[p]<<4) | (cCt[j]<<8) | (cCt[k]<<12) | (cCt[n]<<16) | (cCt[y]<<20);
  ct[2] = cGt[h] | (cGt[p]<<4) | (cGt[j]<<8) | (cGt[k]<<12) | (cGt[n]<<16) | (cGt[y]<<20);
  ct[3] = cTt[h] | (cTt[p]<<4) | (cTt[j]<<8) | (cTt[k]<<12) | (cTt[n]<<16) | (cTt[y]<<20);
  ct[4] = 0;
  ct[5] = 0;

  for (; sc < sl; sc++) {
    p = sc[17];
    ct[0] = (ct[0] << 4) | cAt[p];
    ct[1] = (ct[1] << 4) | cCt[p];
    ct[2] = (ct[2] << 4) | cGt[p];
    ct[3] = (ct[3] << 4) | cTt[p];
    cm = (ct[sc[4]] >> 16) + (ct[sc[3]] >> 12) + (ct[sc[2]] >> 8) +
      (ct[sc[1]] >> 4) + ct[ * sc];

    /* 7 base cloop */

    cv = (cm & 0xf0);
    athresh = 12;
    nch = 0;

    /* exclude the following cloops */
    /* RRnnnNN, NRnnnYN */
    /* NRnnnNN with cstem < 3 Watson-Crick basepairs or equivalent */
    /* RYnnnYN */
    /* NYnnnNN with cstem < 1 Watson-Crick basepair or equivalent */
    /* NYnnnNN with cstem < 2 Watson-Crick basepairs or equivalent */
    /* unless cloop = CTnnnAN */

    if (RI[sc[6]]) {
      if (RI[sc[5]]) goto CLOOP6;
      if (YI[sc[10]]) goto CLOOP6;
      if (cv < 0x60) goto CLOOP6;
    } else {
      if (YI[sc[10]])
        if (RI[sc[5]]) goto CLOOP6;
      if (cv < 0x40) {
        if (cv < 0x20) goto CLOOP6;
        if (sc[5] != Cytosine) goto CLOOP6;
        if (sc[6] != Thymine) goto CLOOP6;
        if (sc[10] != Adenine) goto CLOOP6;
        athresh = 11;
      } else if (cv < 0x70) {
        athresh = 11;
        k = cYI[sc[5]] + cTI[sc[6]] + cRI[sc[10]] + cAI[sc[11]];
        if (sc[6] == Cytosine)
          if (sc[5] == Cytosine)
            k += 16;
          else
        if (sc[5] == Thymine)
          if (sc[11] == Adenine)
            k += 16;
        if (cv == 0x40) {
          if (k < 40) goto CLOOP6;
        } else
        if (cv == 0x50) {
          if (k < 28) goto CLOOP6;
        } else {
          if (k < 20) goto CLOOP6;
          athresh = 9;
        }
      } else
        athresh = (cv < 10) ? 9 : 8;
    }
    chit[0].pos = sc;
    chit[0].stem = 5;
    chit[0].loop = 7;
    chit[0].looppos = sc + 5;
    chit[0].arm = 17;
    chit[0].end = sc + 17;
    chit[0].anticodon = (sc[7] << 4) + (sc[8] << 2) + sc[9];
    if (bp[sc[-1]][sc[17]]) {
      chit[1].pos = sc - 1;
      chit[1].stem = 6;
      chit[1].loop = 7;
      chit[1].looppos = sc + 5;
      chit[1].arm = 19;
      chit[1].end = sc + 18;
      chit[1].anticodon = chit[0].anticodon;
      nch = 2;
    } else nch = 1;

    /* 6 base cloop */
    /* exclude cstem < 4 Watson-Crick basepairs or equivalent */
    /* exclude cloop = RRnnNN */
    /* exclude cloop = NNnnYY */

    CLOOP6:
      if (cloop7) goto CLOOPE;
    if ((cm & 0xf00) >= 0x800) {
      if (!YI[sc[6]])
        if (!YI[sc[5]])
          goto CLOOP8;
      if (!RI[sc[9]])
        if (!RI[sc[10]])
          goto CLOOP8;
      se = sc + 20;
      sg = sc;
      while (sg < se) {
        sf = sg + 5;
        while (sf < (sg + 11)) {
          if ( * sf == * sg)
            if (sf[1] == sg[1])
              if (sf[2] == sg[2])
                if (sf[3] == sg[3])
                  if (sf[4] == sg[4]) {
                    sb = sg + 5;
                    s = sf + 5;
                    i = 0;
                    while (sb < sf)
                      if ( * sb++ != * s++)
                        if (++i > 1) goto NXSEG6;
                    goto CLOOPE;
                  }
          NXSEG6:
            sf++;
        }
        sg++;
      }
      chit[nch].pos = sc;
      chit[nch].stem = 5;
      chit[nch].loop = 6;
      chit[nch].looppos = sc + 5;
      chit[nch].arm = 16;
      chit[nch].end = sc + 16;
      chit[nch++].anticodon = 0;
      if (athresh > 10) athresh = 10;
      if (bp[sc[-1]][sc[16]]) {
        chit[nch].pos = sc - 1;
        chit[nch].stem = 6;
        chit[nch].loop = 6;
        chit[nch].looppos = sc + 5;
        chit[nch].arm = 18;
        chit[nch].end = sc + 17;
        chit[nch++].anticodon = 0;
      }
    }

    /* 8 base cloop */
    /* exclude cstem < 4 Watson-Crick basepairs or equivalent */
    /* exclude cloop = RRnnnnNN */
    /* exclude cloop = NNnnnnYY */

    CLOOP8:
      if ((cm & 0xf) >= 0x8) {
        if (!YI[sc[5]])
          if (!YI[sc[6]])
            goto CLOOPE;
        if (!RI[sc[12]])
          if (!RI[sc[11]])
            goto CLOOPE;
        se = sc + 20;
        sg = sc;
        while (sg < se) {
          sf = sg + 5;
          while (sf < (sg + 11)) {
            if ( * sf == * sg)
              if (sf[1] == sg[1])
                if (sf[2] == sg[2])
                  if (sf[3] == sg[3])
                    if (sf[4] == sg[4]) {
                      sb = sg + 5;
                      s = sf + 5;
                      i = 0;
                      while (sb < sf)
                        if ( * sb++ != * s++)
                          if (++i > 1) goto NXSEG8;
                      goto CLOOPE;
                    }
            NXSEG8:
              sf++;
          }
          sg++;
        }
        chit[nch].pos = sc;
        chit[nch].stem = 5;
        chit[nch].loop = 8;
        chit[nch].looppos = sc + 5;
        chit[nch].arm = 18;
        chit[nch].end = sc + 18;
        chit[nch++].anticodon = 0;
        if (athresh > 10) athresh = 10;
        if (bp[sc[-1]][sc[18]]) {
          chit[nch].pos = sc - 1;
          chit[nch].stem = 6;
          chit[nch].loop = 8;
          chit[nch].looppos = sc + 5;
          chit[nch].arm = 20;
          chit[nch].end = sc + 19;
          chit[nch++].anticodon = 0;
        }
      }

    /* calculate carm energy */

    CLOOPE:
      if (nch < 1) continue;
    for (nc = 0; nc < nch; nc++) {
      s1 = chit[nc].pos;
      cstem = chit[nc].stem;
      cloop = chit[nc].loop;
      s4 = s1 + cstem;
      s2 = s4 + cloop;
      energy = (cloop == 7) ? 0.0 : -4.0;
      energy += cY[ * s4] + cT[s4[1]] + cR[s2[-2]] + cA[s2[-1]];
      if (s4[1] == Cytosine)
        if ( * s4 == Cytosine)
          energy += 2.6;
        else
      if ( * s4 == Thymine)
        if (s2[-1] == Adenine)
          energy += 2.6;
      s2 += cstem;
      stem_energy = bem[ * s1][ * --s2];
      k = neighbour_map[ * s1][ * s2];
      stem_energy += neighbour_em[k][s1[1]][s2[-1]];
      bondtype = btmap[ * s1][ * s2];
      if (bp[ * s1][ * s2]) {
        if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
        else stem_energy += send_em[ * s2][ * s1];
      } else {
        if (assymst[ * s2][ * s1]) stem_energy += mtTERMSTAB;
        else stem_energy += send_em[s2[-1]][s1[1]];
      }
      while (++s1 < s4) {
        if (!wcbp[ * s1][ * --s2]) {
          if (!wcbp[s1[-1]][s2[1]]) {
            for (j = 0; j < mtNTM; j++)
              if ( * s1 == tandemid[j][1])
                if ( * s2 == tandemid[j][3])
                  if (s1[-1] == tandemid[j][0])
                    if (s2[1] == tandemid[j][2]) {
                      stem_energy += tandem_em[j];
                      break;
                    }
            if (s1 < (s4 - 1))
              if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB;
          }
          k = neighbour_map[ * s1][ * s2];
          stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
            neighbour_em[k][s1[1]][s2[-1]]);
        }
        bondtype += btmap[ * s1][ * s2];
        stem_energy += bem[ * s1][ * s2];
      }
      if (!bp[ * --s1][ * s2]) {
        s1--;
        s2++;
      }
      if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
      else stem_energy += send_em[ * s1][ * s2];
      cgcc = bondtype & 0xf;
      if (cgcc <= 0) {
        catc = (bondtype & 0xf0) >> 4;
        if (catc < cstem) energy -= mtGCPENALTY;
      }
      if (cstem == 6) energy += 1.0;
      chit[nc].bondtype = bondtype;
      chit[nc].stem_energy = stem_energy;
      chit[nc].energy = energy + stem_energy;
    }

    /* find tarms */

    nth = 0;
    slm = sc + 61;
    sle = sc + 57;
    sb = sc + 21;
    sg = sc + 16;
    sge = sg + 30;
    slb = sg + 32;
    tem[0] = At[ * slm];
    tem[1] = Ct[ * slm];
    tem[2] = Gt[ * slm];
    tem[3] = Tt[ * slm];
    while (--slm > sle) {
      tem[0] = (tem[0] << 4) | At[ * slm];
      tem[1] = (tem[1] << 4) | Ct[ * slm];
      tem[2] = (tem[2] << 4) | Gt[ * slm];
      tem[3] = (tem[3] << 4) | Tt[ * slm];
    }
    while (slm >= sb) {
      tem[0] = ((tem[0] << 4) | At[ * slm]) & 0xfffff;
      tem[1] = ((tem[1] << 4) | Ct[ * slm]) & 0xfffff;
      tem[2] = ((tem[2] << 4) | Gt[ * slm]) & 0xfffff;
      tem[3] = ((tem[3] << 4) | Tt[ * slm]) & 0xfffff;
      sf = slm + 3;
      if (sf > sge) sf = sge;
      apos2 = slm + 5;
      si = sg;
      s = si + 4;
      r = tem[ * si];
      while (++si < s) r = (r >> 4) + tem[ * si];
      while (si <= sf) {
        if (si < slm)
          r = (r >> 4) + tem[ * si++];
        else {
          si++;
          r = r >> 4;
        }
        q = r & 0xf;
        if (slm > slb) {
          if (q < 5) continue;
          tloop = (int)(slm - si);
        } else {
          if (q < 2) continue;
          if (q < 3) {
            if (!wcbp[si[-5]][apos2[-1]]) continue;
            if (!wcbp[si[-4]][apos2[-2]]) continue;
            tloop = (int)(slm - si);
            if (tloop > 5) continue;
          } else {
            tloop = (int)(slm - si);
            if (q < 4)
              if (!bp[si[-4]][apos2[-2]])
                if (!bp[si[-2]][apos2[-4]]) {
                  if (tloop < 4) continue;
                  if (si[-1] != Guanine) continue;
                  if ( * si != Thymine) continue;
                  if (si[1] != Thymine) continue;
                }
          }
        }
        if (tloop < 7) {
          if (tloop < 2)
            if (tloop <= 0) {
              if (tloop <= -2) {
                if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6;
              } else
              if (bp[si[-3]][apos2[-3]]) {
                tstem = 3;
                tloop += 4;
              } else {
                if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6;
              }
            }
          else {
            if (bp[si[-2]][apos2[-4]]) {
              tstem = 4;
              tloop += 2;
            } else
            if (bp[si[-3]][apos2[-3]]) {
              tstem = 3;
              tloop += 4;
            } else {
              if (!wcbp[si[-5]][apos2[-1]]) continue;
              if (!wcbp[si[-4]][apos2[-2]]) continue;
              tstem = 2;
              tloop += 6;
            }
          } else {
            if (bp[si[-1]][apos2[-5]]) {
              if (q != 4) tstem = 5;
              else {
                if (bp[si[-2]][apos2[-4]]) tstem = 5;
                else {
                  k = GI[si[-3]] + TI[si[-2]] + TI[si[-1]] + CI[ * si];
                  if (k >= 2) {
                    tstem = 3;
                    tloop += 4;
                  } else tstem = 5;
                }
              }
            } else {
              if (bp[si[-2]][apos2[-4]]) {
                tstem = 4;
                tloop += 2;
              } else
              if (bp[si[-3]][apos2[-3]]) {
                tstem = 3;
                tloop += 4;
              } else {
                if (!wcbp[si[-5]][apos2[-1]]) continue;
                if (!wcbp[si[-4]][apos2[-2]]) continue;
                tstem = 2;
                tloop += 6;
              }
            }
          }
          if (tloop < 3)
            if (tstem > 3) {
              tstem--;
              tloop += 2;
            }
        } else {
          if (!bp[si[-1]][apos2[-5]])
            if (!bp[si[-2]][apos2[-4]]) {
              tstem = 3;
              tloop += 4;
            }
          else {
            tstem = 4;
            tloop += 2;
          } else tstem = 5;
        }
        if (tloop > 17)
          if (tstem < 5)
            continue;

        /* calculate tarm energy */

        s1 = si - 5;
        tpos = s1;
        s4 = s1 + tstem;
        s2 = apos2;
        if (tt[ * s1][ * --s2]) {
          energy = mtTSTTSTAB;
          if (tt[ * ++s1][ * --s2]) {
            energy += mtTSTTSTAB;
            bondtype = btmap[ * s1++][ * s2--];
          } else bondtype = 0;
        } else {
          energy = 0.0;
          bondtype = 0;
        }

        /* calculate tstem energy */

        stem_energy = bem[ * s1][ * s2];
        k = neighbour_map[ * s1][ * s2];
        stem_energy += neighbour_em[k][s1[1]][s2[-1]];
        bondtype += btmap[ * s1][ * s2];
        while (++s1 < s4) {
          if (!wcbp[ * s1][ * --s2]) {
            if (!wcbp[s1[-1]][s2[1]]) {
              for (j = 0; j < mtNTM; j++)
                if ( * s1 == tandemid[j][1])
                  if ( * s2 == tandemid[j][3])
                    if (s1[-1] == tandemid[j][0])
                      if (s2[1] == tandemid[j][2]) {
                        stem_energy += tandem_em[j];
                        break;
                      }
              if (s1 < (s4 - 1))
                if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB;
            }
            k = neighbour_map[ * s1][ * s2];
            stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
              neighbour_em[k][s1[1]][s2[-1]]);
          }
          bondtype += btmap[ * s1][ * s2];
          stem_energy += bem[ * s1][ * s2];
        }
        s1--;
        if (tloop < 4) stem_energy += ssend_em[ * s1][ * s2];
        else
        if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
        else stem_energy += send_em[ * s1][ * s2];

        /* compile possible tarms */

        energy += (stem_energy - mtBONDSTAB * (double)(5 - tstem));
        if (energy >= tarmthresh) {
          thit[nth].pos = tpos;
          s1 = tpos + tstem;
          s2 = apos2 - tstem;
          thit[nth].energy = energy - loop_stab[tloop] +
            tG[s1[-1]] + tT[ * s1] + tT[s1[1]] + tC[s1[2]];
          thit[nth].stem_energy = stem_energy;
          thit[nth].bondtype = bondtype;
          thit[nth].stem = tstem;
          thit[nth].loop = tloop;
          thit[nth].end = tpos + 2 * tstem + tloop;
          if (++nth >= mtNTH) {
            fprintf(stderr, "Too many mt-tstem hits\n");
            break;
          }
          if (tstem > 2)
            if (tloop < 10)
              if (gt[s1[-1]][ * s2]) {
                thit[nth].pos = tpos;
                thit[nth].energy = energy - mtBONDSTAB - mtGTBOND -
                  loop_stab[tloop + 2] +
                  tG[s1[-2]] + tT[s1[-1]] + tT[ * s1] + tC[s1[1]];
                thit[nth].stem_energy = stem_energy - mtGTBOND;
                thit[nth].bondtype = bondtype - 0x100;
                thit[nth].stem = tstem - 1;
                thit[nth].loop = tloop + 2;
                thit[nth].end = thit[nth - 1].end;
                if (++nth >= mtNTH) {
                  fprintf(stderr, "Too many mt-tstem hits\n");
                  break;
                }
                if (tstem > 3)
                  if (tloop < 8)
                    if (gt[s1[-2]][s2[1]]) {
                      thit[nth].pos = tpos;
                      thit[nth].energy = energy - 2.0 * mtBONDSTAB - 2.0 * mtGTBOND -
                        loop_stab[tloop + 4] + tG[s1[-3]] +
                        tT[s1[-2]] + tT[s1[-1]] + tC[ * s1];
                      thit[nth].stem_energy = stem_energy - 2.0 * mtGTBOND;
                      thit[nth].bondtype = bondtype - 0x200;
                      thit[nth].stem = tstem - 2;
                      thit[nth].loop = tloop + 4;
                      thit[nth].end = thit[nth - 1].end;
                      if (++nth >= mtNTH) {
                        fprintf(stderr, "Too many mt-tstem hits\n");
                        break;
                      }
                    }
              }
          if (tstem < 5) {
            if (tloop < 11) continue;
            if (tloop > 16) continue;
            if (!wcbp[s1[1]][s2[-2]]) continue;
            bondtype += btmap[ * s1][s2[-1]] + btmap[s1[1]][s2[-2]];
            tstem += 2;
            tloop -= 4;
          } else {
            if (tloop < 9) continue;
            if (wcbp[ * s1][s2[-1]]) {
              if (tloop > 14) continue;
              tstem++;
              tloop -= 2;
              bondtype += btmap[ * s1][s2[-1]];
            } else {
              if (tloop < 11) continue;
              if (tloop > 16) continue;
              if (!wcbp[s1[1]][s2[-2]]) continue;
              bondtype += btmap[ * s1][s2[-1]] + btmap[s1[1]][s2[-2]];
              tstem += 2;
              tloop -= 4;
            }
          }
          thit[nth].pos = tpos;
          s1 = tpos + tstem;
          thit[nth].energy = energy - loop_stab[tloop] +
            tG[s1[-1]] + tT[ * s1] + tT[s1[1]] + tC[s1[2]];
          thit[nth].stem_energy = stem_energy;
          thit[nth].bondtype = bondtype;
          thit[nth].stem = tstem;
          thit[nth].loop = tloop;
          thit[nth].end = thit[nth - 1].end;
          if (++nth >= mtNTH) {
            fprintf(stderr, "Too many mt-tstem hits\n");
            break;
          }
          if (tloop < 9) continue;
          if (!wcbp[ * s1][apos2[-tstem - 1]]) continue;
          if (++tstem > 7) continue;
          if (tloop > 14) continue;
          tloop -= 2;
          thit[nth].pos = tpos;
          s1 = tpos + tstem;
          thit[nth].energy = energy - loop_stab[tloop] +
            tG[s1[-1]] + tT[ * s1] + tT[s1[1]] + tC[s1[2]];
          thit[nth].stem_energy = stem_energy;
          thit[nth].bondtype = bondtype;
          thit[nth].stem = tstem;
          thit[nth].loop = tloop;
          thit[nth].end = thit[nth - 1].end;
          if (++nth >= mtNTH) {
            fprintf(stderr, "Too many mt-tstem hits\n");
            break;
          }
        }
      }
      slm--;
    }

    /* find darms */

    ndh = 0;
    sle = sc - 4;
    slb = sc - 8;
    slm = sc - 1;
    tem[0] = dAt[ * slm];
    tem[1] = dCt[ * slm];
    tem[2] = dGt[ * slm];
    tem[3] = dTt[ * slm];
    while (--slm > sle) {
      tem[0] = (tem[0] << 4) | dAt[ * slm];
      tem[1] = (tem[1] << 4) | dCt[ * slm];
      tem[2] = (tem[2] << 4) | dGt[ * slm];
      tem[3] = (tem[3] << 4) | dTt[ * slm];
    }
    slm1 = slm;
    while (slm > slb) {
      tem[0] = ((tem[0] << 4) | dAt[ * slm]) & 0xffff;
      tem[1] = ((tem[1] << 4) | dCt[ * slm]) & 0xffff;
      tem[2] = ((tem[2] << 4) | dGt[ * slm]) & 0xffff;
      tem[3] = ((tem[3] << 4) | dTt[ * slm]) & 0xffff;
      slm--;
      si = slm - 18;
      s = si + 3;
      r = tem[ * si];
      while (++si < s) r = (r >> 4) + tem[ * si];
      while (si <= slm1) {
        if (si < slm) r = (r >> 4) + tem[ * si++];
        else {
          r = r >> 4;
          si++;
        }
        if ((q = (r & 0xf)) < 6) {
          q += (unsigned int)(TI[si[-6]] + RI[si[-5]]);
          if (q < 6) continue;
        }

        /* calculate darm energy */

        s1 = si - 4;
        dhit[ndh].pos = s1;
        energy = dT[s1[-2]] + dA[s1[-1]];
        dloop = (int)(slm1 - si);
        if (dloop > 2)
          if (bp[si[-1]][ * slm1]) {
            dstem = 4;
            goto EC;
          }
        if (dloop > 0)
          if ((ggstembp[si[-2]][slm[2]]) || (gabp[si[-1]][ * slm1])) {
            dstem = 3;
            dloop += 2;
            energy += mtNOBOND;
            goto EC;
          }
        if (!wcbp[si[-3]][slm[3]]) continue;
        if (!gc[si[-4]][slm[4]]) continue;
        dstem = 2;
        dloop += 4;
        if (dloop > 5) energy += mtNOBOND;
        energy += mtNOBOND;

        EC:
          s2 = slm + 4;
        s4 = s1 + dstem;
        if (!wcbp[s1[1]][s2[-1]])
          if (stemterm[s1[1]][s2[-1]]) energy -= 1.0;
          else
        if (bp[s1[1]][s2[-1]]) energy -= 1.5;
        else energy -= 2.0;

        /* calculate dstem energy */

        stem_energy = bem[ * s1][ * s2];
        k = neighbour_map[ * s1][ * s2];
        stem_energy += neighbour_em[k][s1[1]][s2[-1]];
        bondtype = btmap[ * s1][ * s2];
        if (bp[ * s1][ * s2]) {
          if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
          else stem_energy += send_em[ * s2][ * s1];
          s1++;
          s2--;
        } else {
          s1++;
          s2--;
          if (assymst[s2[1]][s1[-1]]) stem_energy += mtTERMSTAB;
          else stem_energy += send_em[ * s2][ * s1];
        }
        stem_energy += bem[ * s1][ * s2];
        k = neighbour_map[ * s1][ * s2];
        stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
          neighbour_em[k][s1[1]][s2[-1]]);
        bondtype += btmap[ * s1][ * s2];
        while (++s1 < s4) {
          if (!wcbp[ * s1][ * --s2]) {
            if (!wcbp[s1[-1]][s2[1]]) {
              for (j = 0; j < mtNTM; j++)
                if ( * s1 == tandemid[j][1])
                  if ( * s2 == tandemid[j][3])
                    if (s1[-1] == tandemid[j][0])
                      if (s2[1] == tandemid[j][2]) {
                        stem_energy += tandem_em[j];
                        break;
                      }
              if (s1 < (s4 - 1))
                if (!bp[s1[1]][s2[-1]]) stem_energy -= mt3MMSTAB;
            }
            k = neighbour_map[ * s1][ * s2];
            stem_energy += (neighbour_em[k][s1[-1]][s2[1]] +
              neighbour_em[k][s1[1]][s2[-1]]);
          }
          bondtype += btmap[ * s1][ * s2];
          stem_energy += bem[ * s1][ * s2];
        }
        if (!bp[ * --s1][ * s2]) {
          s1--;
          s2++;
        }
        if (dloop < 4) stem_energy += ssend_em[ * s1][ * s2];
        else
        if (assymst[s1[1]][s2[-1]]) stem_energy += mtTERMSTAB;
        else stem_energy += send_em[ * s1][ * s2];

        /* compile possible darms */

        energy += stem_energy;
        dhit[ndh].energy = energy;
        dhit[ndh].stem_energy = stem_energy;
        dhit[ndh].bondtype = bondtype;
        dhit[ndh].stem = dstem;
        dhit[ndh].loop = dloop;
        if (++ndh >= mtND) {
          fprintf(stderr, "Too many mt-dstem hits\n");
          break;
        }
        if (dstem == 4) {
          if (dloop >= 6)
            if (bondtype < 0x1000) {
              s1 = si - 5;
              s2 = slm + 5;
              if (bp[ * s1][ * s2]) {
                dhit[ndh].pos = s1;
                e = 0.5 + bem[ * s1][ * s2];
                dhit[ndh].energy = energy + e;
                if (wcbp[ * s1][ * s2])
                  dhit[ndh].energy += (dT[s1[-2]] + dA[s1[-1]] -
                    dT[s1[-1]] - dA[ * s1]);
                dhit[ndh].stem_energy = stem_energy + e;
                dhit[ndh].bondtype = bondtype + btmap[ * s1][ * s2];
                dhit[ndh].stem = 5;
                dhit[ndh].loop = dloop;
                if (++ndh >= mtND) {
                  fprintf(stderr, "Too many mt-dstem hits\n");
                  break;
                }
              }
            }
        } else
        if (dloop >= 6) {
          s1 = si - 1;
          s2 = slm1;
          if (stemterm[ * s1][ * s2]) {
            dhit[ndh].pos = si - 4;
            dhit[ndh].energy = energy;
            dhit[ndh].stem_energy = stem_energy;
            dhit[ndh].bondtype = bondtype;
            dhit[ndh].stem = 4;
            dhit[ndh].loop = dloop - 2;
            if (++ndh >= mtND) {
              fprintf(stderr, "Too many mt-dstem hits\n");
              break;
            }
          }
        }
        if (dloop >= 4) continue;
        s1 = si - 4 + dstem - 1;
        s2 = s1 + dloop + 1;
        if (bp[ * s1][ * s2]) continue;
        dhit[ndh].pos = si - 4;
        dhit[ndh].energy = energy + 0.001;
        dhit[ndh].stem_energy = stem_energy;
        dhit[ndh].bondtype = bondtype;
        dhit[ndh].stem = dstem - 1;
        dhit[ndh].loop = dloop + 2;
        if (++ndh >= mtND) {
          fprintf(stderr, "Too many mt-dstem hits\n");
          break;
        }
      }
      slm1--;
    }

    /* build darm exclusion map */
    /* 5' astems further from carm than */
    /* mt_DRLmaxlength must match a darm */

    for (i = 3; i <= 30; i++) dposmap[i] = 0;
    sf = sc - mt_DRLmaxlength - 1;
    sld = sf;
    if (ndh > 0) {
      s = dhit[0].pos;
      for (nd = 0; nd < ndh; nd++) {
        se = dhit[nd].pos;
        if (se < s) s = se;
        i = (int)(sc - se);
        if (dposmap[++i] < 1) dposmap[i] = 1;
        dposmap[++i] = 2;
        if (dposmap[++i] < 1) dposmap[i] = 1;
      }
      s -= 4;
      if (s < sf) sf = s;
    }

    /* build tarm exclusion map */
    /* 3' astems further from carm than */
    /* mt_TVRLmaxlength must match a tarm */

    for (i = 17; i <= 62; i++) tendmap[i] = 0;
    s2 = sc + mt_TVRLmaxlength + 17;
    sle = s2;
    if (nth > 0) {
      s = thit[0].end;
      for (nt = 0; nt < nth; nt++) {
        se = thit[nt].end;
        if (se > s) s = se;
        i = (int)(se - sc);
        bondtype = thit[nt].bondtype;
        if (tendmap[i]) {
          if (bondtype < tendmap[i]) tendmap[i] = bondtype;
        } else tendmap[i] = bondtype;
      }
      if (s > s2) s2 = s;
    }

    /* find astems in 3 categories: */
    /* high energy astems close to carm */
    /* high energy astems matching a high energy tarm far from carm */
    /* low energy astem matching a darm and tarm */

    nah = 0;
    sa = sc - 3;
    sg = sf - 6;
    sb = sc + 17;
    se = s2 + 6;
    tem[0] = aAt[ * se];
    tem[1] = aCt[ * se];
    tem[2] = aGt[ * se];
    tem[3] = aTt[ * se];
    while (--se > s2) {
      tem[0] = (tem[0] << 4) | aAt[ * se];
      tem[1] = (tem[1] << 4) | aCt[ * se];
      tem[2] = (tem[2] << 4) | aGt[ * se];
      tem[3] = (tem[3] << 4) | aTt[ * se];
    }
    ti = (int)(se - sc);
    while (se >= sb) {
      tem[0] = ((tem[0] << 4) | aAt[ * se]) & 0xfffffff;
      tem[1] = ((tem[1] << 4) | aCt[ * se]) & 0xfffffff;
      tem[2] = ((tem[2] << 4) | aGt[ * se]) & 0xfffffff;
      tem[3] = ((tem[3] << 4) | aTt[ * se]) & 0xfffffff;
      if (tendmap[ti]) {
        nti = (tendmap[ti] < 0x2000) ? 1 : 0;
      } else {
        if (se > sle) goto ANX;
        nti = -1;
      }
      si = sg;
      r = tem[ * si];
      while (++si < sf) r = (r >> 4) + tem[ * si];
      di = (int)(sc - si);
      while (si < sa) {
        r = (r >> 4) + tem[ * si++];
        if (dposmap[--di]) {
          if (nti <= 0) {
            if (nti < 0)
              if (dposmap[di] < 2) continue;
            if ((av = (r & 0xf)) < athresh) continue;
          }
        } else {
          if (si < sld) continue;
          if (nti < 0) continue;
          if ((av = (r & 0xf)) < athresh) continue;
        }
        if (nah >= mtNA) {
          fprintf(stderr, "Too many mt-astem hits\n");
          break;
        }

        /* predict astem length and calculate astem energy */

        s1 = si - 7;
        s2 = se + 6;
        if (bp[ * s1][ * s2]) {
          astem = 7;
          energy = 0.0;
          ahit[nah].pos1 = s1;
          ahit[nah].pos2 = se;
        } else
        if (ggstemterm[ * s1][ * s2]) {
          astem = 7;
          ahit[nah].pos1 = s1;
          ahit[nah].pos2 = se;
          energy = bem[ * s1++][ * s2--];
        } else {
          energy = bem[ * s1++][ * s2--];
          if (bp[ * s1][ * s2]) {
            astem = 6;
            ahit[nah].pos1 = s1;
            ahit[nah].pos2 = se;
          } else
          if (ggstemterm[ * s1][ * s2]) {
            astem = 6;
            ahit[nah].pos1 = s1;
            ahit[nah].pos2 = se;
            energy += bem[ * s1++][ * s2--];
          } else {
            astem = 5;
            energy += bem[ * s1++][ * s2--];
            ahit[nah].pos1 = s1;
            ahit[nah].pos2 = se;
          }
        }
        ahit[nah].stem = astem;
        bondtype = btmap[ * s1][ * s2];
        energy += bem[ * s1][ * s2];
        k = neighbour_map[ * s1][ * s2];
        energy += neighbour_em[k][s1[1]][s2[-1]];
        energy += bem[ * ++s1][ * --s2];
        k = neighbour_map[ * s1][ * s2];
        energy += (neighbour_em[k][s1[-1]][s2[1]] +
          neighbour_em[k][s1[1]][s2[-1]]);
        bondtype += btmap[ * s1][ * s2];
        while (++s1 < si) {
          if (!wcbp[ * s1][ * --s2]) {
            if (!wcbp[s1[-1]][s2[1]]) {
              for (j = 0; j < mtNTM; j++)
                if ( * s1 == tandemid[j][1])
                  if ( * s2 == tandemid[j][3])
                    if (s1[-1] == tandemid[j][0])
                      if (s2[1] == tandemid[j][2]) {
                        energy += tandem_em[j];
                        break;
                      }
              if (s1 < (si - 1))
                if (!bp[s1[1]][s2[-1]]) energy -= mt3MMSTAB;
            }
            k = neighbour_map[ * s1][ * s2];
            energy += (neighbour_em[k][s1[-1]][s2[1]] +
              neighbour_em[k][s1[1]][s2[-1]]);
          }
          bondtype += btmap[ * s1][ * s2];
          energy += bem[ * s1][ * s2];
        }
        if (!bp[ * --s1][ * s2])
          if (!bp[ * --s1][ * ++s2])
            if (!bp[ * --s1][ * ++s2])
              if (!bp[ * --s1][ * ++s2])
                goto NOST;
        if (assymst[s1[1]][s2[-1]]) energy += mtTERMSTAB;
        NOST:
          ahit[nah].energy = energy;
        ahit[nah].bondtype = bondtype;
        nah++;
      }
      ANX:
        se--;
      ti--;
    }
    if (nah <= 0) continue;

    /* build mttrna genes */
    /* cycle through astems first so that */
    /* GC content is only calculated once per astem */

    thresh = -INACTIVE;
    te.ps = NULL;
    for (na = 0; na < nah; na++) {
      apos2 = ahit[na].pos2;
      apos1 = ahit[na].pos1;
      astem = ahit[na].stem;
      aend1 = apos1 + astem;
      astem8 = (astem == 7) ? (wcbp[apos1[-1]][apos2[7]]) : 0;
      asteme = 0;
      ea = ahit[na].energy;
      abondtype = ahit[na].bondtype;
      agcat = ((abondtype >> 4) + abondtype) & 0xf;

      /* GC content */

      s = apos1;
      aend2 = apos2 + astem;
      nbase = (int)(aend2 - apos1) + 1;
      igc = 0;
      while (s <= aend2) {
        k = * s++;
        if (k >= Cytosine)
          if (k <= Guanine) igc++;
      }
      gcv = 10.0 * (double) igc / (double) nbase;
      if (gcv < 1.0) {
        if (gcv < 0.55) continue;
        ea -= 0.5;
      }
      if (nbase > 60) {
        if (gcv > 6.0) ea -= 2.0 * (gcv - 6.0);
      } else {
        if (gcv > 5.0) ea -= 2.0 * (gcv - 5.0);
      }
      if (gcv > 6.6) {
        ea -= 6.0;
        if (gcv > 7.0) ea -= 6.0;
      }

      /* findout if inside a coding sequence */

      incds = 0;
      i = -1;
      while (++i < ncdsh)
        if (apos1 > cdshit[i].pos1)
          if (aend2 <= cdshit[i].pos2) {
            incds = 1;
            ea -= 2.0;
            break;
          }

      /* cycle through carms that fall between astem */

      nc = -1;
      while (++nc < nch) {
        cpos = chit[nc].pos;
        dloop = (int)(cpos - aend1);
        if (dloop < 3) continue;
        if (dloop > 26) continue;
        cend = chit[nc].end;
        tloop = (int)(apos2 - cend);
        if (tloop < 5) continue;
        cloop = chit[nc].loop;
        cstem = chit[nc].stem;
        clooppos = chit[nc].looppos;
        cloopend = clooppos + cloop;
        carm = chit[nc].arm;
        anticodon = chit[nc].anticodon;
        cbondtype = chit[nc].bondtype;
        acbondtype = abondtype + cbondtype;
        cgcat = ((cbondtype >> 4) + cbondtype) & 0xf;
        ec = ea + chit[nc].energy;

        /* astem,cstem stability (GC bond count) */

        if ((abondtype & 0xf) <= 0)
          if ((cbondtype & 0xf) <= 0) {
            ec -= mtGCPENALTYD;
            if (((cbondtype & 0xf0) >> 4) >= 5) ec += 0.5;
          }

        /* anticodon to astem discriminator base match */

        if (cloop == 7) {
          if (!mt_discrim[par_discrim][anticodon][apos2[astem]])
            if (astem8)
              if (! mt_discrim[par_discrim][anticodon][apos2[8]])
                ec -= 3.0;
          else
          if (astem <= 6) {
            if (!mt_discrim[par_discrim][anticodon][apos2[7]])
              if (astem == 5) {
                if (!mt_discrim[par_discrim][anticodon][apos2[6]]) ec -= 3.0;
              }
            else
              ec -= 3.0;
          } else
            ec -= 3.0;
        }

        /* build TV-replacement loop mttrna genes */

        if (tloop <= mt_TVRLmaxlength) {
          if (!par_tvloop) goto TVN;

          /* astem termination */
          /* (only need to calculate once per astem) */

          if (!asteme) {
            asteme = 1;
            s = aend1 - 1;
            se = apos2;
            while (!bp[ * s][ * se]) {
              if (--s <= apos1) {
                eas = 0.0;
                goto NOST2;
              }
              se++;
            }
            if (!aastemterm[s[1]][se[-1]]) eas = -0.5;
            else {
              eas = 0.0;
              while (se >= apos2) {
                s++;
                se--;
                if (aastemterm[ * s][ * se]) eas += 1.0;
              }
            }
          }

          /* choose darm */

          NOST2:
            energy = 94.0 + ec + eas;
          nd = -1;
          ndi = -1;
          ed = -INACTIVE;
          while (++nd < ndh) {
            dpos = dhit[nd].pos;
            spacer1 = (int)(dpos - aend1);
            if (spacer1 != 2) continue;
            dl = dhit[nd].loop;
            dstem = dhit[nd].stem;
            if (dstem > 4) continue;
            darm = 2 * dstem + dl;
            spacer2 = (int)(cpos - dpos) - darm;

            /* astem,darm,cstem interspacing */

            if (spacer2 < 1) continue;
            e = dhit[nd].energy;
            if (spacer2 > 1) {
              if (spacer2 > 2) continue;
              if (!stembp[ * cpos][cend[-1]]) continue;
              if (tloop > 12) e -= 2.0;
              if ((dhit[nd].bondtype & 0xf) < 1)
                if ((agcat + cgcat + 1) < (cstem + astem)) e -= 3.0;
            } else
            if (dl > 11) {
              if (!RI[cpos[-1]]) e -= 2.0;
            } else {
              if (cpos[-1] == Cytosine) e -= 2.0;
            }

            /* small,large dloop, dstem R motif */

            if (dl < 3) e -= 2.0;
            if (dl > 12) e -= 2.0;
            if (!RI[ * dpos]) e -= 1.0;

            /* darm,tloop tertiary interaction */

            k = 0;
            di = ((dl >= 12) ? 3 : ((dl >= 9) ? 2 : 1));
            tl = (tloop >= 14) ? 5 : ((dl >= 9) ? ((tloop >= 10) ? 4 : 3) : 3);
            if (!ggstackbp[dpos[dstem + di]][cend[tl]]) {
              if (tl > 3) {
                if (!ggstackbp[dpos[dstem + di]][cend[tl - 1]]) e -= 1.5;
                else k++;
              } else
              if (di > 1) {
                if (!ggstackbp[dpos[dstem + di - 1]][cend[tl]]) e -= 1.5;
                else k++;
              } else
                e -= 1.5;
            } else k++;
            if (stemterm[dpos[dstem - 1]][dpos[darm - dstem]]) {
              e -= 0.5;
              if (cend[2] == dpos[dstem - 2]) {
                if (bp[cend[2]][dpos[darm - dstem + 1]]) k++;
              } else {
                if (cend[2] == dpos[darm - dstem + 1])
                  if (bp[cend[2]][dpos[dstem - 2]]) k++;
              }
            } else {
              if (cend[2] == dpos[dstem - 1]) {
                if (!bp[cend[2]][dpos[darm - dstem]]) e -= 0.5;
                else k++;
              } else {
                if (cend[2] != dpos[darm - dstem]) e -= 0.5;
                else
                if (!bp[cend[2]][dpos[dstem - 1]]) e -= 0.5;
                else k++;
              }
            }
            if (cend[1] == * dpos) {
              if (!stackbp[cend[1]][dpos[darm - 1]]) e -= 0.5;
              else k++;
            } else {
              if (cend[1] != dpos[darm - 1]) e -= 0.5;
              else
              if (!bp[cend[1]][ * dpos]) e -= 0.5;
              else k++;
            }

            /* darm stability */

            dstemmotif = wcbp[dpos[1]][dpos[darm - 2]];
            if (spacer2 == 2)
              if ((k < 3) || (dhit[nd].bondtype > 0x200) || (!dstemmotif)) {
                if (abondtype >= 0x10000) e -= 2.0;
                if (dstem > 3) e -= 1.0;
                e -= 0.5;
              }

            /* darm tertiary interactions */

            j = 0;
            b8 = dpos[-2];
            b9 = dpos[-1];
            if (!bp[b8][dpos[dstem]]) e -= 1.0;
            else if (wcbp[b8][dpos[dstem]]) j++;
            if (!bp[b8][dpos[darm - dstem - 1]]) e -= 1.0;
            else if (wcbp[b8][dpos[darm - dstem - 1]]) j++;
            if (!wcbp[dpos[2]][dpos[darm - 3]]) {
              if (!gastembp[b8][dpos[dstem]]) e -= 2.0;
              else if (!gastembp[b8][dpos[darm - dstem - 1]]) e -= 2.0;
              if (!ggstembp[dpos[2]][dpos[darm - 3]]) e -= 1.0;
            } else j++;
            if (!bp[b9][dpos[2]]) {
              if (!bp[b9][dpos[darm - 3]]) e -= 1.0;
              else j++;
            } else j++;

            /* more extensive tertiary interaction between darm,tloop */

            if (dstemmotif) {
              if (k >= 3)
                if (bp[dpos[2]][dpos[darm - 3]]) {
                  if (b8 != Thymine) e += 0.5;
                  if (dl > 3)
                    if (bp[dpos[dstem + 2]][cend[tl + 1]]) e += 0.7;
                    else
                  if (gabp[dpos[dstem + 2]][cend[tl + 1]]) e += 0.5;
                  if (tloop >= 6)
                    if (spacer2 < 2)
                      if (dl >= 3) {
                        di = (dl > 11) ? 2 : 1;
                        if (bp[dpos[dstem + di]][cend[tl]]) {
                          if (chit[nc].stem_energy > -4.8) e += 0.5;
                          if (wcbp[dpos[dstem + di]][cend[tl]])
                            if (gcv > 1.2)
                              if (clooppos[1] == Thymine)
                                if (cbondtype < 0x200)
                                  if ((cbondtype & 0xf) > 0)
                                    if (abondtype < 0x2000) {
                                      e += 1.5;
                                      if (dl > 3)
                                        if (wcbp[dpos[dstem + di + 1]][cend[tl + 1]])
                                          e += 1.0;
                                    }
                        }
                      }
                }
              if (j >= 4) e += 0.25;
            }
            if (e > ed) {
              ed = e;
              ndi = nd;
              ti = k;
            }
          }
          if (ndi < 0) goto TVN;
          energy += ed;
          dpos = dhit[ndi].pos;
          dstem = dhit[ndi].stem;
          dl = dhit[ndi].loop;
          darm = 2 * dstem + dl;
          dbondtype = dhit[ndi].bondtype;
          spacer2 = (int)(cpos - dpos) - darm;
          spacer1 = (int)(dpos - aend1);
          b8 = * aend1;
          b9 = aend1[1];

          /* false positive suppression */

          if (dloop < 15) energy -= 2.0;
          if (cstem > 5) energy -= 1.0;
          if (tloop < 6) energy -= 1.0;
          if (tloop > 12) {
            energy -= 1.0;
            if (agcat < 6) energy -= 2.0;
            if (tloop > 15) energy -= 2.5;
          }
          if (!stackbp[ * dpos][dpos[darm - 1]]) energy -= 1.0;
          if (dstem < 4)
            if (gcv > 1.2)
              if ((dbondtype & 0xf0f) == 0) energy -= 1.5;
          if (b8 != Thymine) {
            if (dl < 4)
              if (abondtype > 0x10000)
                energy -= 1.5;
            if (b8 == Adenine)
              if (YI[cloopend[-2]])
                energy -= 1.0;
          }
          if (dl > 10) {
            if (tloop < 7) energy -= 2.0;
            if (spacer2 > 1) energy -= 2.0;
            if (dhit[ndi].stem_energy < -3.4) energy -= 2.0;
          }
          if (gcv < 2.0)
            if (dbondtype > 0x10000) energy -= 2.0;
          if ((cbondtype & 0xf) < 1)
            if (abondtype > 0x100) {
              if (cgcat < 4)
                energy -= 1.5;
              if (!wcbp[dpos[2]][dpos[darm - 3]]) energy -= 1.0;
            }
          if (b8 != Thymine)
            if ((clooppos[1] != Thymine) ||
              ( * clooppos != Cytosine))
              if (dl > 3)
                if (dbondtype > 0x10000)
                  energy -= 1.0;
          if (!RI[cend[1]])
            if (b9 != Guanine) energy -= 1.0;
            else energy -= 0.5;
          if (b9 == Guanine) {
            if (!RI[ * cend]) energy -= 1.0;
            if (spacer2 != 1) energy -= 3.0;
            else {
              tl = (tloop >= 14) ? 5 : ((dl >= 9) ? ((tloop >= 7) ? 4 : 3) : 3);
              s = dpos + dstem;
              if (!wcbp[s[1]][cend[tl]]) {
                energy -= 2.5;
                if (dl >= 5)
                  if (chit[nc].energy > 2.0)
                    if (wcbp[s[2]][cend[tl]])
                      if (wcbp[s[3]][cend[tl + 1]])
                        energy += 6.0;
              } else
              if (b8 == Thymine)
                if (dl >= 5)
                  if (chit[nc].energy > 2.0)
                    if (wcbp[s[2]][cend[tl + 1]])
                      energy += 3.5;
            }
          } else if (b9 != Adenine) energy -= 3.0;
          if (b8 != Thymine)
            if (b8 == Guanine) {
              if (!RI[dpos[dstem]]) energy -= 1.0;
              else
              if (RI[dpos[darm - dstem - 1]]) energy += 2.0;
            }
          else energy -= 1.0;

          /* carm termination */

          if (assymst[cend[-1]][ * cpos]) energy += 1.0;

          /* CTnnnAA cloop motif */

          energy += CX7[ * clooppos] + AX7[cloopend[-2]];
          if (clooppos[1] == Cytosine) energy -= 2.0;

          /* NNnnnAA cloop motif */

          if (cloopend[-2] == Adenine)
            if (cloopend[-1] == Adenine)
              if (spacer1 == 2)
                if (dbondtype < 0x1000) {
                  if (abondtype < 0x100) energy += 1.0;
                  else
                  if (cbondtype < 0x100) energy += 1.0;
                }

          /* global stem damage level */

          bondtype = acbondtype + dbondtype;
          i = (int)((bondtype >> 16) & 0xf);
          j = (int)((bondtype >> 12) & 0xf);
          k = (int)((bondtype >> 8) & 0xf);
          if (k > 0)
            if (i > 0) {
              k += (i + j);
              if (k > 5) energy -= 1.0 * (double)(k - 5);
            }

          /* global stem stability (GC bond count) */

          gcc = bondtype & 0xf;
          if (gcc < 2) {
            if (ti >= 2) {
              if (cbondtype < 0x100)
                if ((cbondtype & 0xf) > 0) goto NGCC1;
              if (ti >= 3)
                if (cgcat >= 4) {
                  if ((cbondtype & 0xf) > 0) goto NGCC1;
                  if (cbondtype < 0x100) goto NGCC2;
                }
            }
            energy -= (double)(3 - gcc);
            NGCC2:
              if (gcc < 1) {
                if (agcat < 5) energy -= 2.0;
                if (bondtype > 0x10000) energy -= 1.5;
              }
          }
          NGCC1:

            /* global stability */
            /* (stem stability,dloop-tloop tertiary interaction,dloop size) */

            if (abondtype > 0x1000)
              if (ti < 3) {
                if (chit[nc].stem_energy < -6.0)
                  energy -= 1.5;
                if (dl > 9)
                  if (((dbondtype + cbondtype) & 0xf) < 1)
                    energy -= 1.0;
              }

          /* tloop,dloop tertiary interaction */
          /* (alternative dloop position) */

          if (bondtype < 0x1000)
            if (b8 == Thymine)
              if (RI[b9])
                if (dl > 4)
                  if (!bp[cend[3]][dpos[dstem + 1]])
                    if (bp[cend[3]][dpos[dstem + 2]])
                      energy += 0.5;

          /* "near perfect" TV-loop mttRNA: */
          /* darm-tloop tertiary interaction,low global stem damage, */
          /* TR motif at b8-9, good astem,darm,carm interspacing */

          if (ti >= 2)
            if (agcat >= 6)
              if (cbondtype < 0x100)
                if (dbondtype < 0x100)
                  if (RI[b9])
                    if (b8 == Thymine)
                      if ((abondtype & 0xf) > 0)
                        if ((dbondtype & 0xf) > 0)
                          if (spacer1 == 2)
                            if (spacer2 == 1)
                              energy += 1.5;

          /* find exceptions */

          if (energy < dthresh) {
            if (!mtxdetect) goto TVN;
            if (incds) goto TVN;
            if (energy < (thresh - 7.0)) goto TVN;
            if (energy < (dthresh - 7.0)) goto TVN;
            if (nbase > 68) goto TVN;
            if (abondtype > 0x20100) goto TVN;
            if (dl > 9) {
              if (dl > 10) goto TVN;
              if (dstem < 4) goto TVN;
              if (dbondtype > 0x100) goto TVN;
            }
            if (dstem > 4) goto TVN;
            if (b9 != Adenine) {
              if (b9 != Guanine) goto TVN;
              if (cbondtype > 0x100) goto TVN;
              if (dbondtype > 0x200) goto TVN;
            }
            if (cloop != 7) goto TVN;
            if (YI[cloopend[-2]]) goto TVN;
            if (b8 == Thymine) {
              if (apos2[-1] == Thymine)
                if (apos2[-2] == Thymine)
                  if (tloop < 8)
                    if (tt[aend1[-1]][ * apos2])
                      if (wcbp[dpos[2]][dpos[darm - 3]])
                        if (((dbondtype + cbondtype) & 0xf) > 0)
                          energy += 3.0;
            } else
            if (b8 == Adenine) {
              if (apos2[-1] == Adenine)
                if (apos2[-2] == Adenine) {
                  if (assymat[aend1[-1]][ * apos2])
                    if (assymat[apos2[1]][aend1[-2]]) energy += 2.0;
                  if (agcat >= 5)
                    if (cgcat >= 4)
                      if (dbondtype < 0x100)
                        if (at[aend1[-1]][ * apos2])
                          if (at[apos2[1]][aend1[-2]])
                            energy += 1.0;
                }
              if (ti >= 3)
                if (cgcat >= 4)
                  if (agcat >= 4)
                    if ((cbondtype & 0xf) > 0)
                      if ((abondtype & 0xf) > 1)
                        if (dbondtype < 0x200)
                          if (wcbp[dpos[1]][dpos[darm - 2]])
                            if (clooppos[1] == Thymine)
                              if (YI[ * clooppos])
                                if (RI[cloopend[-2]])
                                  if (RI[cloopend[-1]])
                                    energy += 5.0;
            }
            if (bondtype < 0x100) {
              if (spacer2 == 1)
                if ( * clooppos == Cytosine)
                  if (clooppos[1] == Thymine)
                    if (cloopend[-2] == Adenine)
                      if (cloopend[-1] == Adenine)
                        energy += 2.0;
            } else {
              if (spacer2 == 1) {
                if (b8 == Thymine)
                  if (dl > 3)
                    if (dbondtype < 0x200) {
                      if (cbondtype < 0x100) {
                        if (!bp[dpos[dstem + 1]][cend[3]])
                          if (bp[dpos[dstem + 1]][cend[4]])
                            energy += 2.0;
                        if (dbondtype < 0x100)
                          if (abondtype < 0x20000)
                            if (ti >= 2)
                              if (dstem >= 3)
                                if (tloop < 13)
                                  if ((cbondtype & 0xf) > 0)
                                    energy += 4.0;
                      }
                    }
                else
                if (dstem > 3)
                  if (dbondtype < 0x300) {
                    if (bondtype < 0x10000)
                      if (ti >= 3)
                        if ((acbondtype & 0xf) > 0)
                          if (wcbp[dpos[2]][dpos[darm - 3]])
                            energy += 4.0;
                  }
                if (tloop < 8) {
                  if (dbondtype < 0x200) {
                    if (cbondtype < 0x100)
                      if (ti >= 2) {
                        if (wcbp[dpos[dstem + 1]][cend[3]]) {
                          if (b8 == Thymine)
                            if (abondtype < 0x3000)
                              energy += 5.0;
                          if (agcat >= 5)
                            if (gcv > 1.2)
                              if (RI[cloopend[-1]])
                                energy += 7.0;
                        }
                        if (dbondtype < 0x100)
                          if (agcat >= 6)
                            if (YI[ * clooppos])
                              if (clooppos[1] == Thymine)
                                if (RI[cloopend[-2]])
                                  if (RI[cloopend[-1]])
                                    energy += 2.0;
                      }
                    if (cbondtype < 0x300)
                      if (ti >= 3)
                        if (abondtype < 0x2000)
                          if ((dbondtype & 0xf) > 0)
                            if ((acbondtype & 0xf) > 0)
                              if (ahit[na].energy >= -7.0)
                                if (dstem >= 4)
                                  energy += 3.0;
                  }
                  if (dbondtype < 0x300)
                    if (cgcat >= 4)
                      if (abondtype < 0x2000)
                        if (ahit[na].energy >= -7.0)
                          if (cbondtype < 0x10000)
                            if ((cbondtype & 0xf) > 0)
                              if (cstem < 6)
                                if (ti >= 3)
                                  energy += 4.0;
                }
              }
              if (tloop > 8)
                if (agcat >= 6)
                  if (cbondtype < 0x100)
                    if ((cbondtype & 0xf) > 0)
                      if (b8 == Thymine)
                        if (wcbp[dpos[dstem + 1]][cend[3]])
                          if (wcbp[dpos[1]][dpos[darm - 2]])
                            energy += 7.0;
            }
            if (dbondtype < 0x100)
              if (cgcat >= 4)
                if (agcat >= 5)
                  if (wcbp[dpos[1]][dpos[darm - 2]])
                    if ((cbondtype & 0xf) > 0)
                      if ((abondtype & 0xf) > 0)
                        if ((dbondtype & 0xf) > 0)
                          energy += 0.5;
            if (cbondtype < 0x100)
              if (dbondtype < 0x200)
                if (agcat >= 5)
                  if (b8 == Thymine)
                    if (tloop < 8)
                      if (wcbp[dpos[1]][dpos[darm - 2]])
                        if (wcbp[dpos[2]][dpos[darm - 3]])
                          if ((cbondtype & 0xf) > 0)
                            if ((abondtype & 0xf) > 0)
                              if ((dbondtype & 0xf) > 0)
                                if (clooppos[1] == Thymine)
                                  if (YI[ * clooppos])
                                    if (RI[cloopend[-2]])
                                      energy += 3.0;
            if (energy < dthresh) goto TVN;
            energy -= (0.9 * (energy - dthresh) + 5.0);
          }

          /* remember fully formed TV-loop replacement mttRNA gene */
          /* if threshold reached */

          if (energy < thresh) goto TVN;
          te.energy = energy;
          thresh = energy;
          te.ps = apos1;
          te.dstem = dstem;
          te.dloop = dl;
          te.spacer1 = spacer1;
          te.spacer2 = spacer2;
          te.cstem = cstem;
          te.cloop = cloop;
          k = astem + spacer1 + darm + spacer2;
          te.anticodon = k + cstem + 2;
          te.nintron = 0;
          te.intron = 0;
          te.var = 0;
          te.varbp = 0;
          te.tstem = 0;
          te.tloop = tloop;
          te.nbase = k + carm + tloop;

          /* build D-replacement loop mttrna genes */

          TVN:
            if (tloop < 10) continue;
        }
        if (dloop > mt_DRLmaxlength) goto DN;
        if (gcv < 1.2) goto DN;
        energy = 91.0 + ec;

        /* CCnnnAA cloop */

        if (clooppos[1] == Cytosine) {
          if ( * clooppos != Cytosine) goto DN;
          if (cloopend[-2] != Adenine) goto DN;
          if (cloopend[-1] != Adenine) goto DN;
          energy -= 1.0;
        }

        /* choose tarm */

        nt = -1;
        nti = -1;
        et = -INACTIVE;
        while (++nt < nth) {
          tl = thit[nt].loop;
          if (tl > 11) continue;
          if (thit[nt].end != apos2) continue;
          tpos = thit[nt].pos;
          tstem = thit[nt].stem;

          /* var loop (3-7 bases long) */

          var = (int)(tpos - cend);
          if (var < 3) continue;
          e = thit[nt].energy;
          if (var > 5) {
            if (var > 7) continue;
            if (tl < 7) continue;
            e -= 1.0;
            if ((dloop < 10) || (tstem < 4))
              e -= 2.0 * (double)(var -5);
          }

          /* tloop RA or RG motif */

          s = tpos + tstem;
          k = 0;
          n = 0;
          i = 0;
          while ((j = tloopa[tl][i++]) >= 0)
            if (s[j] == Adenine) {
              k = 1;
              if (dloop >= 3)
                if (tl > 3) {
                  b57 = s[j - 1];
                  if (RI[b57] || (tl < 5)) {
                    if (bp[b57][aend1[0]]) {
                      e += 1.5;
                      n = 1;
                      break;
                    }
                    if (bp[b57][aend1[1]]) {
                      e += 1.5;
                      n = 1;
                      break;
                    }
                    if (dloop > 10)
                      if (bp[b57][aend1[2]]) {
                        e += 1.5;
                        n = 1;
                        break;
                      }
                  }
                }
            }
          if (!k) {
            i = 0;
            while ((j = tloopa[tl][i++]) >= 0)
              if (s[j] == Guanine)
                if (RI[s[j - 1]]) {
                  k = 1;
                  break;
                }
            if (j < 0) e -= ((tl > 5) ? 2.0 : 1.0);
          }

          /* tertiary interaction between tloop and start of dloop */

          ti = (tl > 5) ? 1 : ((dloop > 5) ? 1 : 0);
          di = (dloop > 5) ? 2 : 1;
          if (stackbp[aend1[di]][s[ti]]) e += 1.0;

          /* tloop GTTC motif */

          i = (s[-1] == Guanine) ? 1 : 0;
          if (tl >= 5) {
            ti = i + TI[ * s] + TI[s[1]] + CI[s[2]];
            if (n)
              if (!i)
                if (TI[ * s])
                  if (TI[s[1]])
                    if (AI[s[2]])
                      if (tl >= 7)
                        ti++;
            if ((i > 0) || (ti >= 3)) e += (double) ti;
          } else {
            ti = i + TI[ * s] + TI[s[1]];
            if ((i > 0) || (ti >= 2)) e += (double) ti;
          }
          if (e > et) {
            et = e;
            nti = nt;
            tc = k;
          }
        }

        if (nti < 0) goto DN;
        energy += et;
        tpos = thit[nti].pos;
        tstem = thit[nti].stem;
        tl = thit[nti].loop;
        tbondtype = thit[nti].bondtype;
        var = (int)(tpos - cend);

        /* tertiary interaction between b48(=tpos[-1]) and dloop */

        b48 = tpos[-1];
        if (dloop <= 7) {
          if (YI[b48]) tc++;
          else energy -= 1.0;
        } else {
          i = 0;
          while ((j = dloopi[dloop][i++]) >= 0)
            if (assymagbp[b48][aend1[j]]) {
              tc++;
              break;
            }
          if (j < 0) energy -= 1.0;
        }

        /* large dloop, large tloop */

        if (dloop > 7) {
          if (tl >= 6)
            if (tc < 2)
              energy -= 2.0;
          if (tstem < 3)
            energy -= 1.0;
        }

        /* carm termination */

        s = cpos - 1;
        se = cend;
        if (cstem > 5) {
          s++;
          se--;
        }
        if (!stackbp[ * s][ * se]) energy -= 1.0;
        se = cpos - 3;
        if (!bp[cend[-1]][ * cpos]) {
          if (assymst[cend[-1]][ * cpos]) {
            if (dloop < 5) se++;
            energy += 1.5;
          } else if (dloop < 13) se++;
        } else {
          if (cstem > 5) {
            if (dloop < 13) se++;
          } else if (dloop < 5) se++;
        }

        /* tertiary interaction between tloop and dloop near carm */

        s = tpos + tstem;
        if (tl >= 5) {
          ti = (tl >= 10) ? 4 : ((tl >= 7) ? 3 : 2);
          b57 = s[ti];
          if (!gabp[ * se][b57]) energy -= 2.0;
          else {
            k = (var > 3) ? 2 : ((var > 1) ? 1 : 0);
            if (bp[cend[k]][b57]) energy += 1.0;
          }
        }

        /* R motif at end of tstem */

        if (!RI[s[-1]]) energy -= 2.0;

        /* large tloop */

        if (tl > 9)
          if (tbondtype > 0x200) energy -= 2.0;

        /* dloop,var,tloop T repeat motif */
        /* present in some nematode D-loop replacement tRNA-Ser genes */

        if (dloop >= 4) {
          k = 1;
          se = aend1;
          while (se < cpos)
            if ( * se++ == Thymine) k++;
          if (k >= dloop) {
            if (var >= 3) {
              se = cend;
              while (se < tpos)
                if ( * se++ == Thymine) k++;
              if (k >= (var +dloop)) {
                energy += 3.0;
                se = s + ((tl > 5) ? 5 : tl);
                while (s < se)
                  if ( * s++ != Thymine) break;
                if (s >= se) energy += 5.5;
              }
            }
          }
        }

        /* astem stability  */

        if (ea < -6.1)
          if (tl > 4) {
            if ( * s == Thymine)
              if (s[-1] == Guanine)
                if (s[1] == Thymine)
                  goto NASI;
            if (ea > -8.3)
              if ( * clooppos == Cytosine)
                if (clooppos[1] == Thymine)
                  if (cloopend[-2] == Adenine)
                    if (cloopend[-1] == Adenine)
                      goto NASI;
            energy -= 3.0;
          }
        NASI:

          /* cstem stability (GC bond count) */

          bondtype = acbondtype + tbondtype;
        if ((cbondtype & 0xf) < 1)
          if ((bondtype & 0xf) < 3) energy -= 1.0;

        /* cloop CTnnnAA motif */

        if (bondtype >= 0x400)
          energy += CX[ * clooppos] + TX[clooppos[1]] +
          AXX[cloopend[-1]] + AXX37[cloopend[-2]];
        else
          energy += CX[ * clooppos] + TX[clooppos[1]] +
          AX[cloopend[-1]] + AX37[cloopend[-2]];

        /* large dloop */

        if (dloop >= 9) {
          k = tloop - dloop - 4;
          if (k < 0)
            if (bondtype >= 0x1000) energy += (double) k;
          if (dloop >= 12) {
            if (dloop >= 14) energy -= 2.0;
            else
            if (tstem < 6) energy -= ((dloop >= 13) ? 2.0 : 1.0);
          }
        }

        /* small dloop, small tarm */

        if (dloop <= 10)
          if (tstem < 3)
            if (ea > -2.6)
              if (tl <= 7)
                if (cgcat >= 4)
                  if (gc[ * tpos][apos2[-1]])
                    if (gc[tpos[1]][apos2[-2]])
                      if (gcv > 1.2)
                        if ((abondtype & 0xf) > 0)
                          if ((cbondtype & 0xf) > 0)
                            energy += (4.5 + (mtBONDSTAB - 0.5) * (double)(5 - tstem));

        /* global stem damage level */

        i = (int)((bondtype >> 16) & 0xf);
        j = (int)((bondtype >> 12) & 0xf) + i;
        k = (int)((bondtype >> 8) & 0xf);
        if (tstem > 3) {
          if ((k > 0) || (tl > 9))
            if ((j > 0) || (k > 5)) {
              n = j + k;
              if ((s[-1] != Guanine) || ( * s != Thymine) ||
                (s[1] != Thymine) || (tstem < 5))
                if (n > 4) energy -= 2.0 * (double)(n - 4);
            }
        } else {
          n = j + k;
          if (n > 3) energy -= 2.0 * (double)(n - 3);
        }

        /* long tstem with tloop GTT motif */

        if (s[-1] == Guanine)
          if ( * s == Thymine)
            if (s[1] == Thymine)
              if (tstem >= 6)
                if (tbondtype < 0x100)
                  energy += 1.5;

        /* find exceptions */

        if (energy < tthresh) {
          if (!mtxdetect) goto DN;
          if (incds) goto DN;
          if (energy < (thresh - 13.5)) goto DN;
          if (energy < (tthresh - 13.5)) goto DN;
          if (k > 1) {
            if (i > 2) goto DN;
            if (k > 4)
              if (i > 1) goto DN;
          }
          if (nbase > 70) goto DN;
          if (var > 4) {
            if (var > 5) goto DN;
            if (var > tl) goto DN;
          }
          if (tstem < 4)
            if ((agcat + cgcat + 2) < (astem + cstem))
              goto DN;
          if (tl > 9) goto DN;
          if (dloop > 13) goto DN;
          if (!YI[ * clooppos]) goto DN;
          if ((abondtype & 0xf) < 2) {
            if ((abondtype & 0xf) < 1) goto DN;
            if (cbondtype > 0x200)
              if (tbondtype > 0x100)
                if (abondtype > 0x200)
                  goto DN;
          }
          if ((tbondtype & 0xf) < 1) {
            if ((acbondtype & 0xf) < 1) goto DN;
            if (acbondtype > 0x200) goto DN;
          }
          if ((dloop + 19) < tloop) goto DN;
          if (gcv > 5.5) goto DN;
          tgcat = ((tbondtype >> 4) + tbondtype) & 0xf;
          if ((tgcat + 2) < tstem) goto DN;
          if (cloop != 7) goto DN;
          if (bp[ * cpos][cend[-1]])
            if (bp[cpos[-1]][ * cend])
              if (bp[cpos[-2]][cend[1]])
                energy += 2.0;
          if (bondtype < 0x20000)
            if (thit[nti].stem_energy > -4.6)
              if (tstem >= 4)
                if ((tstem >= 5) || (s[-1] == Guanine))
                  if (stackbp[cpos[1]][cend[-2]])
                    if (stackbp[ * cpos][cend[-1]])
                      if (stackbp[cpos[-1]][ * cend]) {
                        energy += 1.5;
                        if (s[-1] == Guanine)
                          if ( * s == Thymine)
                            if (s[1] == Thymine)
                              energy += 1.0;
                        if (agcat >= 6) energy += 0.5;
                      }
          if (tc > 0)
            if (tstem >= 5)
              if (var < 6)
                if (var > 2) {
                  if (acbondtype < 0x100) energy += 5.0;
                  else
                  if ((abondtype + tbondtype) < 0x100)
                    energy += 3.0;
                  else
                  if (cloopend[-2] == Thymine)
                    if (cloopend[-1] == Thymine)
                      if (dloop > 7)
                        if (tbondtype < 0x100)
                          if (!tt[ * tpos][apos2[-1]])
                            if ((agcat + cgcat) >= 10)
                              energy += 13.5;
                }
          if (s[-1] == Guanine)
            if ( * s == Thymine)
              if (s[1] == Thymine)
                if ((tstem >= 5) || (s[2] == Cytosine)) {
                  energy += 1.5;
                  if (tstem >= 5)
                    if (tbondtype < 0x1000)
                      if (s[2] == Cytosine) {
                        if (abondtype < 0x10000) {
                          if ( * clooppos == Cytosine)
                            if (clooppos[1] == Thymine)
                              if (cloopend[-2] == Adenine)
                                if (cloopend[-1] == Adenine)
                                  energy += 3.0;
                          if (tbondtype < 0x200)
                            if (bondtype < 0x10000)
                              if (tl == 7)
                                if (s[4] == Adenine)
                                  energy += 4.0;
                        }
                      }
                  else
                  if (tbondtype < 0x200)
                    if ((tbondtype & 0xf) >= 2)
                      if ( * clooppos == Cytosine)
                        if (clooppos[1] == Thymine)
                          if (cloopend[-2] == Adenine)
                            if (cloopend[-1] == Adenine)
                              energy += 1.0;
                }
          if (tstem >= 4)
            if (tbondtype < 0x100)
              if (cbondtype < 0x200)
                if (agcat >= 5) energy += 1.5;
          if (energy > tthresh) energy = tthresh;
          if (ea > -1.8) energy += 3.0;
          else
          if (abondtype < 0x60) energy += 1.5;
          else
          if (acbondtype < 0x200) energy += 0.75;
          if ( * clooppos == Cytosine)
            if (cloopend[-2] == Adenine)
              if (cloopend[-1] == Adenine) {
                if (tstem >= 5)
                  if (tbondtype < 0x100)
                    if (clooppos[1] == Thymine) {
                      energy += 3.0;
                      if (tstem >= 6) energy += 1.0;
                    }
                else
                if (clooppos[1] == Cytosine)
                  energy += 1.0;
                if (tc >= 2)
                  if (clooppos[1] == Thymine)
                    if (bondtype < 0x1000)
                      if (tstem >= 4)
                        if (var < 6)
                          if (var > 2)
                            energy += 3.0;
              }
          if (cbondtype < 0x100)
            if (agcat >= 5)
              if (tc > 0)
                if (clooppos[1] == Thymine)
                  if (YI[ * clooppos])
                    if (RI[cloopend[-2]])
                      if (RI[cloopend[-1]])
                        if (tbondtype < 0x100)
                          energy += 4.0;
                        else
          if (agcat >= 6)
            if ((tgcat + 1) >= tstem)
              if (tstem >= 4)
                energy += 4.0;
          if (bondtype < 0x1000) {
            energy += 0.5;
            if (bondtype < 0x200) energy += 0.75;
          }
          if (energy < tthresh) goto DN;
          energy -= (3.0 + 0.9 * (energy - tthresh));
        }

        /* mammalian cloop motif constraint */

        if (par_discrim == MAMMAL_MT) {
          s1 = clooppos;
          s2 = s1 + cloop;
          r = * s1++;
          while (s1 < s2) r = (r << 4) + * s1++;
          if (r != clmotif[0])
            if (r != clmotif[1])
              if (r != clmotif[2])
                energy -= 5.0;
        }

        /* remember fully formed D-loop replacement mttRNA gene */
        /* if threshold reached */

        if (energy < thresh) goto DN;
        te.energy = energy;
        thresh = energy;
        te.ps = apos1;
        te.spacer1 = 0;
        te.spacer2 = 0;
        te.dstem = 0;
        te.dloop = dloop;
        te.cstem = cstem;
        te.cloop = cloop;
        te.anticodon = astem + dloop + cstem + 2;
        te.nintron = 0;
        te.intron = 0;
        te.var =
          var;
        te.varbp = 0;
        te.tstem = tstem;
        te.tloop = tl;
        te.nbase = astem + dloop + carm +
          var +
          2 * tstem + tl;

        /* build fully formed cloverleaf mttRNA genes */

        DN:
          if (dloop < 10) continue;

        /* choose tarm */

        nt = -1;
        nti = -1;
        et = -INACTIVE;
        while (++nt < nth) {
          tend = thit[nt].end;
          if (tend != apos2) continue;
          e = thit[nt].energy;
          tpos = thit[nt].pos;
          tstem = thit[nt].stem;

          /* GT motif on tloop */

          s = tpos + tstem;
          if ( * s == Thymine)
            if (s[-1] == Guanine)
              if (tstem >= 5)
                if (!stackbp[ * tpos][tend[-1]]) {
                  e += 0.5;
                  if (!bp[tpos[1]][tend[-2]]) e += 0.5;
                }

          /* large var loop */

          var = (int)(tpos - cend);
          if (var > 5) {
            ev = (double)(var -5);
            if (tstem < 5) e -= 3.0 * ev;
            else e -= (0.5 + 0.5 * ev);

            /* allow large var loop if tarm looks nuclear */
            /* (GTTC motif, very large var loop base-pairing) */

            if (var > 9) {
              if ((thit[nt].bondtype & 0xf) < 1) e -= 1.0;
              e -= (0.25 * (double)(var -8));
              if ( * s == Thymine)
                if (s[-1] == Guanine)
                  if (s[1] == Thymine)
                    if (s[2] == Cytosine)
                      e += 4.0;
              if (var > 17) {
                if (var > 25) continue;
                e += 0.5 * vloop_stability(cend,
                  var, & varbp);
              }
            }
          }

          /* small var loop */

          if (var < 3) {
            if (tstem > 5)
              if (s[-1] != Guanine)
                e -= 0.5;
            if (var < 2) {
              if (var < 1) {
                if (var < 0) continue;
                if (tstem < 4)
                  if (thit[nt].stem_energy < -4.0)
                    continue;
              }
              e -= 3.0;
            }
          }
          if (e > et) {
            et = e;
            nti = nt;
          }
        }

        if (nti < 0) continue;
        tpos = thit[nti].pos;
        tstem = thit[nti].stem;
        tl = thit[nti].loop;
        tarm = 2 * tstem + tl;
        var = (int)(tpos - cend);
        b48 = tpos[-1];
        tbondtype = thit[nti].bondtype;
        bondtype = acbondtype + tbondtype;
        ti = (int)(((bondtype >> 16) & 0xf) + ((bondtype >> 12) & 0xf) +
          ((bondtype >> 8) & 0xf));

        /* choose darm */

        nd = -1;
        ndi = -1;
        ed = -INACTIVE;
        while (++nd < ndh) {
          dl = dhit[nd].loop;
          dstem = dhit[nd].stem;
          darm = 2 * dstem + dl;
          dpos = dhit[nd].pos;
          e = dhit[nd].energy;

          /* spacing between astem,darm,carm */

          spacer1 = (int)(dpos - aend1);
          spacer2 = (int)(cpos - dpos) - darm;
          if (spacer1 < 2) {
            if (spacer1 < 1) continue;
            if (dstem < 3) continue;
            if (dl > 12) e -= 2.0;
            if (astem < 7) e -= 1.0;
            if (spacer2 != 2) {
              if (spacer2 < 1) continue;
              if (spacer2 > 2) continue;
              if ((abondtype & 0xf) < 1)
                if ((dhit[nd].bondtype & 0xf) < 1)
                  e -= 0.5;
              if (var > 7) e -= 1.0;
              if (dl > 12) e -= 1.0;
              if (cloop != 7) e -= 2.0;
              if (cstem < 6) e -= 3.6;
              else e -= 0.5;
            } else {
              if (cstem > 5) continue;
              s = cpos;
              se = cend - 1;
              while (!bp[ * s][ * se]) {
                s++;
                se--;
              }
              if (!stemterm[s[-1]][se[1]]) e -= 0.5;
              e -= 0.8;
            }
          } else {
            if (spacer1 > 2) {
              if (spacer1 > 3) continue;
              if (dstem > 4) continue;
              if (dstem < 3) continue;
              if (tl > 15) continue;
              if (astem < 7) e -= 1.0;
              if (ti > 4) e -= 1.0;
              if (cloop != 7) e -= 2.0;
              if (tbondtype > 0x2000)
                if (!RI[tpos[tstem - 1]]) e -= 2.0;
              e -= 1.0;
              if (spacer2 != 1) e -= 0.5;
              else
              if (dhit[nd].bondtype < 0x100)
                if (var >= 3)
                  if (var <= 5)
                    if (tstem >= 3) {
                      e += 1.0;
                      if (agcat >= 5)
                        if (wcbp[ * aend1][ * apos2])
                          if (!bp[aend1[-1]][ * apos2])
                            if (bp[b48][dpos[dstem + 1]])
                              e += 0.5;
                    }
            }
            if (spacer2 > 1) {
              if (spacer2 > 2) continue;
              if (astem < 7)
                if (spacer1 == 2)
                  e -= 1.0;
              if (cloop != 7) e -= 2.0;
              if (ea < -5.8) e -= 2.0;
              e -= 2.5;
              if (bp[b48][dpos[dstem + 1]]) {
                if (dhit[nd].bondtype < 0x1000)
                  if (wcbp[dpos[1]][dpos[darm - 2]])
                    if (wcbp[dpos[2]][dpos[darm - 3]])
                      if (var < 6)
                        if (dl > 3)
                          e += 2.0;
              } else e -= 1.0;
            } else
            if (spacer2 < 1) {
              if (spacer2 < 0) continue;
              if (var > 6) continue;
              if (dstem > 4) continue;
              if (dhit[nd].stem_energy < -4.3) continue;
              if (astem < 7)
                if (spacer1 == 2)
                  e -= 1.0;
              if (cloop != 7) e -= 2.0;
              e -= mtBONDSTAB;
            }
            if (cstem > 5)
              if ((!gt[ * cpos][cend[-1]]) || astem8) e -= mtBONDSTAB;
          }

          /* very large or very small dloop */

          if (dl < 3) e -= 2.0;
          if (dl > 11) {
            if (dl > 14) e -= 2.0;
            else
            if (dl > 13) {
              if (dhit[nd].bondtype >= 0x100) e -= 2.0;
              else e -= 1.0;
            } else
            if (dl > 12) {
              if (dhit[nd].bondtype >= 0x1000) e -= 2.0;
              else e -= 1.0;
            } else
            if (dhit[nd].bondtype >= 0x10000) e -= 2.0;
          }

          /* tertiary interactions in darm */

          b8 = dpos[-2];
          b9 = dpos[-1];
          if (dl > 2) {
            if (dl > 5)
              if (!stackbp[dpos[dstem + 1]][b48]) e -= 1.0;
            if (!stackbp[b8][dpos[dstem]]) e -= 0.25;
            if (!stackbp[b8][dpos[dstem + dl - 1]]) e -= 0.25;
          }
          if (!bp[b9][dpos[2]])
            if (!bp[b9][dpos[darm - 3]])
              e -= 1.0;

          /* TR motif at b8-9 */

          if (RI[b9]) {
            if (b8 == Thymine)
              if (spacer1 == 2)
                if (ti < 6)
                  if (((bondtype & 0xf) > 2) || (bondtype < 0x1000) ||
                    ((tbondtype < 0x100) && (tstem > 3)))
                    if ((cbondtype & 0xf) < 5)
                      if (stembp[dpos[1]][dpos[darm - 2]])
                        if (var < 6)
                          if (var > 2) e += 1.5;
                          else
            if (tstem > 3)
              if (cloopend[-2] == Adenine)
                e += 1.5;
          } else {
            e -= 1.0;
            if (b9 == Thymine)
              if (spacer1 == 2) e -= 2.0;
          }
          if (e > ed) {
            ed = e;
            ndi = nd;
          }
        }

        if (ndi < 0) continue;
        energy = 100.0 + ec + ed + et;
        dl = dhit[ndi].loop;
        dstem = dhit[ndi].stem;
        darm = 2 * dstem + dl;
        dpos = dhit[ndi].pos;
        dbondtype = dhit[ndi].bondtype;
        spacer1 = (int)(dpos - aend1);
        spacer2 = (int)(cpos - dpos) - darm;
        b8 = dpos[-2];

        /* tertiary structure interaction between tloop and dloop */

        if (tl >= 3)
          if (dl >= 4) {
            di = (dl < 7) ? (darm - dstem - 2) : (darm - dstem - 3);
            ti = (tl < 9) ? (tstem + 2) : ((tl < 13) ? (tstem + 3) : (tstem + 5));
            if (ggbp[dpos[di]][tpos[ti]])
              if (ggbp[dpos[di - 1]][tpos[ti - 1]]) {
                energy += 2.0;
                if (spacer1 != 2)
                  if (spacer2 != 2)
                    if (dstem < 4)
                      if (tl > 7)
                        if (bp[dpos[di + 1]][tpos[ti + 1]])
                          energy += 4.0;
                if (ea > -2.5)
                  if (wcbp[dpos[1]][dpos[darm - 2]])
                    if (wcbp[dpos[2]][dpos[darm - 3]])
                      energy += 3.0;
              }
            if (tl > 10)
              if (dl > 10)
                energy -= 1.0;
          }
        else
        if (dl == 3)
          if (wcbp[dpos[dstem + 1]][b48]) energy += 1.0;

        /* small darm and tarm */

        if (tloop <= 18)
          if (tarm <= 13)
            if (dl <= 8)
              if (spacer1 == 2)
                if (spacer2 == 1)
                  if (abondtype < 0x1000)
                    if (tbondtype < 0x100)
                      if (dbondtype < 0x200) {
                        et = (mtBONDSTAB - 0.5) * (double)(5 - tstem) +
                          0.1 * (double)(7 - tl);
                        ed = mtBONDSTAB * (double)(4 - dstem);
                        energy += (0.8 * (et + ed));
                      }

        /* GTTC motif on tloop */

        s = tpos + tstem;
        if (tl < 5)
          if (tl < 2) energy += G[s[-1]];
          else {
            et = (G[s[-1]] + T[ * s] + T[s[1]]);
            if (tl > 3)
              if (bp[ * s][s[tl - 1]]) {
                e = (G[ * s] + T[s[1]] + T[s[2]]);
                if (e > et) et = e;
              }
            if (tstem < 5) {
              e = (G[s[-2]] + T[s[-1]] + T[ * s] + C[s[1]]);
              if (e > et) et = e;
            }
            energy += et;
          }
        else energy += (G[s[-1]] + T[ * s] + T[s[1]] + C[s[2]]);

        /* long astem */

        if (astem8)
          if (bp[apos1[0]][apos2[6]])
            if (bp[apos1[1]][apos2[5]])
              if (bp[apos1[2]][apos2[4]])
                if (bp[apos1[3]][apos2[3]])
                  energy += hbem[apos1[-1]][apos2[7]];

        /* false positive supression */

        if (!RI[cend[0]]) energy -= 1.0;
        if (!RI[cpos[-1]]) energy -= 1.0;
        if (tarm < (var +3)) energy -= 2.0;
        if (gcv < 1.5)
          if (dbondtype > 0x10000) energy -= 2.0;
        if (tarm > 27) {
          energy -= 1.0;
          if (spacer2 != 1) energy -= 1.0;
        }
        if (dstem < 3) {
          if (var > 5) energy -= 1.0;
          if (tloop > (dloop + 8)) energy -= 0.5;
        }
        if (b8 != Thymine)
          if (dl > 3)
            if (dbondtype > 0x100)
              if ((b8 == Cytosine) ||
                (dbondtype > 0x10000))
                if ( * clooppos != Cytosine)
                  if (!wcbp[dpos[dstem + 1]][b48])
                    energy -= 1.0;

        /* high GC false positive suppression */

        if (gcv >= 5.1) {
          if ((abondtype & 0xf) >= 4) {
            s1 = apos1;
            s2 = apos2 + astem;
            n = 0;
            while (--s2 >= apos2)
              if (gc[ * s1++][ * s2]) {
                if (++n >= 4) {
                  energy -= 2.0;
                  break;
                }
              }
            else n = 0;
          }
          if ((dbondtype & 0xf) >= 4) energy -= 3.0;
          if ((cbondtype & 0xf) >= 5) energy -= 3.5;
          if ((tbondtype & 0xf) >= tstem) energy -= 4.0;
        }

        /* global stem damage level */

        tc = tstem + dstem;
        dtbondtype = dbondtype + tbondtype;
        mabondtype = dtbondtype + cbondtype;
        bondtype = acbondtype + dtbondtype;
        if (bondtype < 0x100) energy += 0.5;
        if ((dtbondtype & 0xf) < 1) {
          energy -= 1.0;
          if (tc >= 10) energy -= 2.0;
          if ((bondtype & 0xf) < 3)
            if (nbase > 75) energy -= 1.0;
        }
        i = (int)((bondtype >> 16) & 0xf);
        j = (int)((bondtype >> 12) & 0xf) + i;
        k = (int)((bondtype >> 8) & 0xf) + j;
        ti = (tc > 6) ? 5 : ((tc > 5) ? 4 : 3);
        if (k > ti) {
          ev = (double)(k - ti);
          energy -= 0.5 * ev;
          if (cbondtype > 0x10000)
            if (tstem < 5)
              energy -= ev;
          if (i > 0)
            if (k > 8)
              energy -= 1.5 * (double)(k - 8);
        }

        /* low GC false positive supression */

        if (gcv < 3.5)
          if ((bondtype & 0xf) < 2) {
            if ((bondtype & 0xf) < 1) energy -= 1.0;
            if (dl > 3)
              if (var > 2)
                if (!wcbp[dpos[dstem + 1]][b48]) energy -= 1.0;
          }

        /* small variable loop */

        if (var < 3) {
          if (dloop > 18) {
            if (dloop > (tloop + 2)) energy -= 1.0;
            if (tloop > 20)
              if ((((dtbondtype >> 4) + dtbondtype) & 0xf) < 6)
                energy -= 2.0;
          }
          if (astem < 7) {
            energy -= 1.0;
            if (agcat >= 5)
              if (bondtype < 0x300)
                if (gcv > 1.2)
                  if (gcv < 5.0)
                    energy += 2.0;
          }
        } else

          /* NNNNNAA cloop */

          if (cloopend[-2] == Adenine)
            if (cloopend[-1] == Adenine)
              if (spacer1 > 1)
                if ((dbondtype < 0x2000) || (dloop > mt_DRLmaxlength)) {
                  if (abondtype < 0x100) energy += 1.0;
                  else
                  if (cbondtype < 0x100) energy += 1.0;
                  else
                  if (tstem >= 5)
                    if (tbondtype < 0x100) {
                      energy += 1.0;
                      if ( * clooppos == Cytosine)
                        if (clooppos[1] == Thymine)
                          if (dbondtype < 0x100)
                            energy += 0.5;

                      if (cgcat >= 3)
                        if ((tbondtype & 0xf) > 0)
                          if (ggbp[dpos[dstem + 1]][b48])
                            if (wcbp[dpos[1]][dpos[darm - 2]])
                              if (tl < 10)
                                if (spacer1 == 2)
                                  if (spacer2 == 1)
                                    if (dl > 2)
                                      if (var >= 2)
                                        if (var < 6) {
                                          if (agcat >= 6) energy += 3.0;
                                          else
                                          if (agcat >= 5)
                                            if (cgcat >= 4)
                                              if (dbondtype < 0x100)
                                                if ( * s == Thymine)
                                                  if (s[-1] == Guanine)
                                                    if (s[1] == Thymine)
                                                      energy += 3.0;
                                        }
                    }
                }

        /* large tloop */

        if (tl > 12) {
          if (tbondtype > 0x10000) energy -= 2.0;
          if (agcat < 5)
            if (spacer1 != 2)
              if (spacer2 != 1)
                energy -= 1.0;
        }

        /* find exceptions */

        if (energy < dtthresh) {
          if (!mtxdetect) continue;
          if (incds) continue;
          if (energy < (thresh - 12.0)) continue;
          if (energy < (dtthresh - 12.0)) continue;
          if (nbase > 75) continue;
          if (dstem > 4) continue;
          if (dstem < 3) continue;
          if (astem < 7)
            if (acbondtype > 0x21000)
              continue;
          if (var > 5) {
            if (var > 6) continue;
            if (tarm < 12) continue;
          }
          if (gcv <= 1.2) {
            if (gcv < 0.9) continue;
            if ((mabondtype & 0xf) < 1) continue;
          }
          if (tl > 9) {
            if (tl > 13) continue;
            if (!wcbp[dpos[1]][dpos[darm - 2]])
              continue;
          }
          if (dl > 7) {
            if (bondtype > 0x20000)
              if (dloop > (tloop + 4)) continue;
            if (dl > 10) {
              if (dl > 12)
                if (abondtype > 0x1000)
                  continue;
              if (tbondtype > 0x200) continue;
              if (tt[ * tpos][apos2[-1]]) continue;
              if (var > 5) continue;
              if (dloop > (tloop + 8))
                if (bondtype > 0x10000) continue;
              if (astem < 7) continue;
            }
          }
          if (RI[clooppos[1]]) continue;
          b9 = dpos[-1];
          if (cstem >= 6) {
            if (cbondtype > 0x200) continue;
            if (var < 3) continue;
            if (YI[b9]) continue;
          }
          if (cloop != 7) continue;
          if (par_discrim == MAMMAL_MT) continue;
          if (mabondtype < 0x400) {
            if ((b8 == Thymine) || (mabondtype < 0x300))
              if (ea < -5.45)
                if (chit[nc].stem_energy > -3.2)
                  if (dbondtype < 0x200)
                    if (spacer1 > 1)
                      if ((spacer2 == 1) || (mabondtype < 0x100))
                        if ((spacer1 < 3) || (tstem > 3) ||
                          (tbondtype < 0x100))
                          if ((spacer1 < 3) ||
                            ((var > 2) && (var < 6) && (tbondtype < 0x2000) &&
                              (tl < 10)))
                            if (dstem < 5)
                              if (var >= 2)
                                if (dl > 2)
                                  if (tl < 15)
                                    if ((b8 != Cytosine) ||
                                      ( * clooppos == Cytosine))
                                      if (RI[b9])
                                        if ( * clooppos != Adenine)
                                          if (clooppos[1] == Thymine)
                                            if (RI[cloopend[-2]]) {
                                              s1 = apos1;
                                              s2 = apos2 + astem;
                                              n = 0;
                                              while (--s2 >= apos2)
                                                if (wcbp[ * s1++][ * s2]) {
                                                  if (++n >= 3) break;
                                                }
                                              else n = 0;
                                              if (n >= 3) {
                                                energy += 3.0;
                                                if ((abondtype & 0xf) > 0) energy += 2.0;
                                                if (bp[dpos[dstem + 1]][b48])
                                                  if (wcbp[dpos[1]][dpos[darm - 2]])
                                                    if (var <= 5)
                                                      energy += 1.0;
                                              }
                                              if (dtbondtype < 0x200)
                                                if (agcat < 2)
                                                  if (wcbp[dpos[dstem + 1]][b48])
                                                    if (wcbp[dpos[1]][dpos[darm - 2]])
                                                      if (wcbp[dpos[2]][dpos[darm - 3]])
                                                        if (gcv > 1.2)
                                                          if (var <= 5)
                                                            if (tstem >= 3)
                                                              if (dstem >= 3)
                                                                if (tl > 3)
                                                                  if (tl < 9)
                                                                    if (dl < 9)
                                                                      if (spacer1 == 2)
                                                                        energy += 10.0;
                                            }
            if ((tbondtype & 0xf) > 0)
              if (mabondtype < 0x300) {
                if (mabondtype < 0x100) {
                  if ((spacer1 < 3) || (tstem > 2))
                    if (var > 0)
                      if (YI[ * clooppos])
                        if ((spacer2 > 0) || (clooppos[1] == Thymine))
                          energy += 2.5;
                } else
                if ((dbondtype & 0xf) > 0)
                  if (b9 != Cytosine)
                    if (var <= 7)
                      if (spacer2 == 1)
                        if (tarm < 22)
                          if (gcv > 1.2)
                            if (dstem >= 4) {
                              if (tstem >= 5)
                                energy += 5.0;
                              else
                              if (tstem >= 3)
                                if (tbondtype < 0x100)
                                  energy += 1.0;
                            }
                else
                if (tstem >= 5)
                  energy += 1.0;
              }
            else
            if ((dbondtype & 0xf) > 0) {
              if (tstem >= 5)
                if (s[-1] == Guanine)
                  if ( * s == Thymine)
                    if (s[1] == Thymine)
                      if ( * clooppos == Cytosine)
                        if (clooppos[1] == Thymine)
                          if (cloopend[-2] == Adenine)
                            if (cloopend[-1] == Adenine)

                              energy += 1.0;
              if (bondtype < 0x1000)
                if (cbondtype < 0x100)
                  energy += 1.0;
            }
          }

          if (tstem >= 5)
            if ( * clooppos == Cytosine) {
              if (dl > 3)
                if (dtbondtype < 0x200)
                  if ((tbondtype & 0xf) > 0)
                    if (clooppos[1] == Thymine) {
                      if (clooppos[2] == Thymine)
                        if (clooppos[3] == Adenine)
                          if (clooppos[4] == Cytosine)
                            if (clooppos[5] == Adenine)
                              if (cloop == 7)
                                energy += 0.5;

                      if (cgcat >= 4)
                        if (wcbp[dpos[1]][dpos[darm - 2]])
                          if (bp[dpos[dstem + 1]][b48])
                            if (tl < 10)
                              if (var < 6)
                                if (spacer1 == 2)
                                  if (spacer2 == 1)
                                    if (dstem >= 3)
                                      energy += 3.0;

                    }

              if (clooppos[1] == Cytosine)
                if (clooppos[2] == Cytosine)
                  if (clooppos[3] == Adenine)
                    if (clooppos[4] == Thymine)
                      if (s[-1] == Guanine)
                        if ( * s == Thymine)
                          if (s[1] == Thymine)
                            energy += 1.0;
            }

          if (RI[b9]) {
            if (b8 == Thymine) {
              if (clooppos[1] == Thymine) {
                if (cloopend[-2] == Adenine) {

                  if (wcbp[dpos[1]][dpos[darm - 2]]) {
                    if ( * clooppos == Cytosine) {
                      if (abondtype < 0x200)
                        energy += 1.0;

                      if (bondtype < 0x10000)
                        if (dtbondtype < 0x200)
                          if (agcat >= 3)
                            if (cgcat >= 4)
                              if (tl < 10)
                                if (var < 6)
                                  if (spacer1 == 2)
                                    if (spacer2 == 1)
                                      if (tstem >= 3)
                                        energy += 3.0;
                    }

                    if (tstem >= 5)
                      if (s[-1] == Guanine)
                        if ( * s == Thymine)
                          if (s[1] == Thymine) {
                            energy += 1.0;
                            if (tl >= 5)
                              if (spacer1 == 2)
                                if (spacer2 == 1)
                                  if (tbondtype < 0x100)
                                    if (wcbp[dpos[dstem + 1]][b48])
                                      energy += 3.0;
                          }

                    if (tstem >= 3)
                      if (tl < 10)
                        if (spacer1 == 2)
                          if (spacer2 == 1)
                            if (RI[cloopend[-1]])
                              if (dl > 2)
                                if (var >= 2)
                                  if (var < 6)
                                    if (ggbp[dpos[dstem + 1]][b48]) {
                                      if (dtbondtype < 0x100) {
                                        energy += 3.5;
                                        if ((bondtype & 0xf00) == 0)
                                          if ( * clooppos == Cytosine)
                                            energy += 1.5;
                                      }

                                      if (bondtype < 0x10000)
                                        if (tstem > 2)
                                          if (tbondtype < 0x200) energy += 2.5;

                                      if (abondtype < 0x100)
                                        if (wcbp[dpos[2]][dpos[darm - 3]])
                                          energy += 3.0;

                                      if (tbondtype < 0x100)
                                        if (agcat >= 6)
                                          if (tstem >= 5)
                                            if ((tbondtype & 0xf) > 0)
                                              if (RI[cloopend[-1]])
                                                if (cgcat >= 4)
                                                  energy += 2.0;
                                    }
                    else
                    if (!ggstembp[ * tpos][apos2[-1]])
                      if (wcbp[dpos[dstem + 1]][ * tpos])
                        energy += 1.5;
                  }

                  if ((abondtype & 0xf) < 1)
                    if (abondtype < 0x100)
                      if (gcv > 1.2)
                        if (dl > 3)
                          if (bp[dpos[dstem + 1]][b48])
                            if (spacer1 == 2)
                              if (spacer2 == 1)
                                if ( * clooppos == Cytosine)
                                  energy += 5.0;

                  if (cbondtype < 0x100)
                    if (tbondtype < 0x100)
                      if (tstem >= 3)
                        if (dl > 3)
                          if (var < 6)
                            if (bp[dpos[dstem + 1]][b48])
                              if (spacer1 == 2)
                                if (spacer2 == 1)
                                  energy += 2.5;

                }

                if (stembp[dpos[dstem + 1]][b48]) {

                  if ( * clooppos == Thymine)
                    if (cloopend[-2] == Guanine)
                      if (clooppos[2] == Guanine)
                        if (clooppos[3] == Thymine)
                          if (clooppos[4] == Guanine)
                            if (dl > 2)
                              energy += 1.0;

                  if (cbondtype < 0x100)
                    if (dbondtype < 0x10000)
                      if (wcbp[dpos[1]][dpos[darm - 2]])
                        if (var < 6)
                          if (tstem >= 3)
                            if (gcv >= 1.2)
                              if (dl > 3)
                                energy += 1.0;

                  if (tstem >= 5)
                    if (dtbondtype < 0x200)
                      if ( * clooppos == Cytosine)
                        if (spacer1 == 2)
                          if (spacer2 == 1)
                            if (RI[cloopend[-2]])
                              energy += 0.5;
                }

              }

              if (tstem > 2)
                if (tarm < 28)
                  if (spacer1 == 2)
                    if (spacer2 == 1)
                      if (dl > 3)
                        if (j < 1)
                          if (k > ti)
                            if (ggstembp[dpos[dstem + 1]][b48])
                              energy += 2.5;

              if (dtbondtype < 0x100)
                if ((tbondtype & 0xf) > 0)
                  if (bp[dpos[dstem + 1]][b48])
                    if (b9 == Adenine)
                      if ((dbondtype & 0xf) > 0) energy += 2.0;
                      else
              if (spacer2 == 1)
                energy += 0.5;
              if (cloopend[-2] == Adenine)
                if (cloopend[-1] == Adenine)
                  if (cbondtype < 0x2000)
                    if (spacer1 > 1)
                      if (dl > 2)
                        if (var < 6)
                          energy += 0.75;

            }

            if (var > 2)
              if (dl > 2) {
                if (cbondtype < 0x200)
                  if (((mabondtype & 0xf) > 3) || (bondtype < 0x1000)) {
                    if (bp[dpos[dstem + 1]][b48])
                      energy += 1.0;

                    if (cbondtype < 0x100)
                      if (dbondtype < 0x100)
                        if (bp[b8][dpos[dstem]])
                          if (bp[b8][dpos[darm - dstem - 1]])
                            if (var < 6)
                              if (tstem >= 3)
                                if (tl < 10)
                                  if (spacer1 == 2)
                                    if (spacer2 == 1)
                                      if (clooppos[1] == Thymine)
                                        if (cloopend[-2] == Adenine)
                                          energy += 3.0;
                  }

                if (clooppos[1] == Thymine)
                  if (RI[cloopend[-2]]) {
                    if ( * clooppos == Cytosine) {

                      if (dtbondtype < 0x200)
                        if (agcat >= 3)
                          if (cgcat >= 4)
                            if (var < 6)
                              if (tstem >= 3)
                                if (tl < 10)
                                  if (spacer1 == 2)
                                    if (spacer2 == 1) {
                                      if (abondtype > 0x20000)
                                        if (bp[dpos[dstem + 1]][b48])
                                          energy += 7.0;
                                      if (agcat >= 6) energy += 2.0;
                                    }

                      if ((bondtype & 0xf00) == 0)
                        if (gcv > 5.0)
                          if (s[-1] == Guanine)
                            if ( * s == Thymine)
                              if (tstem >= 5)
                                if (var < 6)
                                  if (tl < 10)
                                    if (spacer1 == 2)
                                      if (spacer2 == 1)
                                        energy += 2.0;

                      if (abondtype < 0x100)
                        if (cbondtype < 0x10000)
                          if (bp[dpos[dstem + 1]][b48])
                            if (cgcat >= 4)
                              if (tstem >= 3)
                                if (var < 6)
                                  if (tstem >= 3)
                                    if (tl < 10)
                                      if (spacer1 == 2)
                                        if (spacer2 == 1)
                                          energy += 1.5;
                    }

                    if (dtbondtype < 0x100)
                      if (agcat >= 4)
                        if (cgcat >= 4)
                          if (var < 6)
                            if (tstem >= 3)
                              if (tl < 10) {
                                if (spacer1 == 2) {
                                  if (abondtype < 0x3000)
                                    if (stackbp[dpos[dstem + 1]][b48])
                                      energy += 3.0;
                                  if (b8 == Thymine)
                                    if (s[-1] == Guanine)
                                      if ( * s == Thymine)
                                        if (s[1] == Thymine)
                                          energy += 3.5;
                                }
                                if (agcat >= 6)
                                  if (YI[ * clooppos])
                                    if (s[-1] == Guanine)
                                      if ( * s == Thymine)
                                        if ((dtbondtype & 0xf) > 0)
                                          energy += 3.0;
                              }

                    if (mabondtype < 0x10000)
                      if (dtbondtype < 0x400)
                        if (agcat >= 5)
                          if (cgcat >= 3)
                            if (tl < 10)
                              if (var < 6)
                                if (spacer1 == 2)
                                  if (spacer2 == 1) {

                                    if (dtbondtype < 0x200)
                                      if (cbondtype < 0x300)
                                        if (bondtype < 0x10000)
                                          if (tstem >= 3)
                                            energy += 1.0;

                                    if (tstem >= 5)
                                      if (s[-1] == Guanine)
                                        energy += 4.0;
                                  }
                  }
              }

          } else
          if (bondtype < 0x10000)
            if (mabondtype < 0x500)
              if (dbondtype < 0x100)
                if (b8 == Thymine)
                  if (agcat >= 4)
                    if (clooppos[1] == Thymine)
                      if (cloopend[-2] == Adenine)
                        if (cloopend[-1] == Adenine)
                          if (spacer1 == 2)
                            if (spacer2 == 1)
                              if (dstem >= 3)
                                if (tstem >= 5)
                                  if (dl > 2)
                                    if (tl < 10)
                                      if (var < 6)
                                        energy += 7.0;

          if (agcat >= 5)
            if (cgcat >= 4)
              if ((acbondtype & 0xf) >= 3) {
                if (tbondtype < 0x100)
                  if ((dbondtype & 0xf) > 0)
                    if ((((dbondtype >> 4) + dbondtype) & 0xf) >= 3)
                      if (wcbp[dpos[dstem + 1]][b48])
                        if (b8 == Thymine)
                          if (RI[b9])
                            if (clooppos[1] == Thymine)
                              if (YI[ * clooppos])
                                if (RI[cloopend[-2]])
                                  if (RI[cloopend[-1]])
                                    if (spacer1 == 2)
                                      if (spacer2 == 1)
                                        if (dl > 2)
                                          if (tl < 10)
                                            if (var < 6)
                                              if (var > 2)
                                                energy += 6.0;
                if (cgcat >= 5)
                  if (abondtype < 0x10000)
                    if (bp[dpos[dstem + 1]][b48])
                      if (clooppos[1] == Thymine)
                        if (YI[ * clooppos])
                          if (RI[cloopend[-2]])
                            if (dl > 2)
                              if (tl < 10)
                                if (var < 6)
                                  if (var > 2)
                                    energy += 6.0;
              }

          if (energy >= dtthresh)
            energy -= (0.9 * (energy - dtthresh) + 5.0);
          else continue;
        }

        /* remember fully formed mttRNA gene if threshold reached */

        if (energy < thresh) continue;
        te.energy = energy;
        thresh = energy;
        te.ps = apos1;
        te.spacer1 = spacer1;
        te.dstem = dstem;
        te.dloop = dl;
        te.spacer2 = spacer2;
        te.cstem = cstem;
        te.cloop = cloop;
        te.var =
          var;
        te.varbp = (var > 17) ? varbp : 0;
        te.tstem = tstem;
        te.tloop = tl;
        k = astem + spacer1 + darm + spacer2;
        te.anticodon = k + cstem + 2;
        te.nintron = 0;
        te.intron = 0;
        te.nbase = k + carm +
          var +2 * tstem + tl;
      }
    }

    /* for highest energy mttRNA gene */
    /* decide astem length, look for NCCA acceptor tail */
    /* and calculate total length */

    if (te.ps) {
      apos2 = te.ps + te.nbase;

      /* store mttRNA gene if there are no */
      /* higher energy overlapping mttRNA genes */

      te.start = (long)(te.ps - seq);

      gs.push_back(make_trna(te));
    }
  }

  gs = best_hit(gs);

  return (gs);
}
