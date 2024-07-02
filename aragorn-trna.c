#include <stdio.h>
#include <stdlib.h>

#define INACTIVE        2.0e+35

#define tRNA    0
#define tmRNA   1
#define rRNA    3
#define CDS     4
#define NS      6    /* should be one more than number of types of gene */

#define TERM   -1

#define Adenine         0
#define Cytosine        1
#define Guanine         2
#define Thymine         3
#define AMBIG           4
#define NOBASE          5

#define ASTEM2_EXT      9
#define MINTSTEM_DIST   (17 + ASTEM2_EXT)
#define MAXTSTEM_DIST   (26 + ASTEM2_EXT)
#define MAXINTRONLEN    3000
#define MINCTRNALEN     62
#define MAXCTRNALEN     110
#define MINTRNALEN      (MINCTRNALEN + 1)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM2_EXT)
#define MAXETRNALEN     (MAXTRNALEN + MAXINTRONLEN)
#define VARMAX          26
#define VARMIN          3
#define VARDIFF         23                /* VARMAX - VARMIN */
#define MINTPTSDIST     50
#define MAXTPTSDIST     321
#define MINTPDIST       50
#define MAXTPDIST       250
#define MAXTAGDIST      102

#define NA          MAXINTRONLEN
#define ND          100
#define NTH         3000
#define NC          5000
#define ATBOND      2.5
#define GCBOND      3.0


typedef struct { int* seq;
                 long size; } dna_sequence;

typedef struct { int *ps;
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
                 int tpe; } gene;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 double energy; } trna_loop;

typedef struct { int *pos;
                 int *end;
                 int stem;
                 int loop;
                 double energy; } trna_dloop;

typedef struct { int *pos1;
                 int *pos2;
                 int stem;
                 double energy; } trna_astem;

typedef struct { int trna;
                 int tmrna;
                 int mtrna;
                 int cloop7;
                 int extastem;
                 int aataildiv;
                 int sp1max;
                 int sp2min;
                 int sp2max;
                 int maxintronlen;
                 int minintronlen;
                 int ifixedpos;
                 int tmstrict;
                 int loffset;
                 int roffset;
                 double threshlevel;
                 double trnathresh;
                 double ttscanthresh;
                 double ttarmthresh;
                 double tdarmthresh;
                 double tastemthresh;
                 double tascanthresh;
                 double mttthresh;
                 double mtdthresh;
                 double mtdtthresh;
               } csw;

typedef struct genes {
  gene *gene;
  struct genes *next;
} genes;

genes* add_gene(genes* gs, gene g){
  genes* new_root = (genes*)malloc(sizeof(genes));
  new_root->next = gs;
  gene* new_gene = (gene*)malloc(sizeof(gene));
  *new_gene = g;
  new_root->gene = new_gene;
  return new_root; 
}


/* Basepair matching matrices */

int bp[6][6] = { { 0,0,0,1,1,0 },
                 { 0,0,1,0,1,0 },
                 { 0,1,0,1,1,0 },
                 { 1,0,1,0,1,0 },
                 { 1,1,1,1,1,0 },
                 { 0,0,0,0,0,0 } };

int vbp[6][6] =
 { { 0,0,1,4,4,0 },
   { 0,0,4,0,4,0 },
   { 1,4,0,2,4,0 },
   { 4,0,2,0,4,0 },
   { 4,4,4,4,4,0 },
   { 0,0,0,0,0,0 } };

double bem[6][6] =
 { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
   { -0.428,-2.144, GCBOND,-2.144, 0.000, 0.000 },
   { -2.144, GCBOND,-2.144, 1.286, 0.000, 0.000 },
   {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
   {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
   {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };


/* LIBRARY */

double vloop_stability(int *sb, int var, int *varbp){
  int e,stem,vstem,loop,*sn,*sen,*pos1,*pos2,*se,*sc,*sd,*sf,*s;
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

// find all tstems in the sequence
//
// tstems are added as trna_loop structs to the hit array
// the hit array is preallocated to some arbitrary size (3000 by default)
//
// modifies hit
// returns the number of tstems found
//
// T-Loop with canonical numbering:
//
//                           (60) (59)
//  (65) (64) (63) (62) (61)          (58)
//    |    |    |    |    |            (57)
//  (49) (50) (51) (52) (53)          (56)
//                           (54) (55)
//
int find_tstems(
  dna_sequence *seq, // sequence
  trna_loop hit[],   // array of trna_loops (currently all NULL)
  int nh,            // hit length
  csw *sw            // state
) {

  int *s = seq->seq;
  int ls = seq->size;

  int i,r,c,tstem,tloop,ithresh1;
  int *s1,*s2,*se,*ss,*si,*sb,*sc,*sf,*sl,*sx,*tem;
  double ec,energy,penalty,thresh2;
  static double bem[6][6] = {
  //    A        C       G       T
     { -2.144,  -0.428, -2.144,  ATBOND, 0.000, 0.000 }, // A
     { -0.428,  -2.144,  GCBOND, -2.144,  0.000, 0.000 }, // C
     { -2.144,   GCBOND, -2.144,  1.286,  0.000, 0.000 }, // G
     {  ATBOND, -2.144,  1.286, -0.428,  0.000, 0.000 }, // T
     {  0.000,   0.000,  0.000,  0.000,  0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000,  0.000, 0.000 } };

  static double A[6] = { 2.0,0.0,0.0,0.0,0.0,0.0 };
  static double C[6] = { 0.0,2.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,2.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,2.0,0.0,0.0 };
  static int tem_trna[6] = { 0x0100, 0x0002, 0x2000, 0x0220, 0x0000, 0x0000 };
  static int tem_tmrna[6] = { 0x0100, 0x0002, 0x2220, 0x0220, 0x0000, 0x0000 };
  i = 0;
  tem = (sw->tmrna || (sw->threshlevel < 1.0))?tem_tmrna:tem_trna;
  ithresh1 = (int)sw->ttscanthresh;
  thresh2 = sw->ttarmthresh;
  ss = s + sw->loffset;
  si = ss + 4 - 1;
  sl = s + ls - sw->roffset + 5 + 3;
  r = tem[*si++];
  r = (r >> 4) + tem[*si++];
  r = (r >> 4) + tem[*si++];

  while (si < sl) {
    r = (r >> 4) + tem[*si++];
    if ((c = (r & 0xF)) < ithresh1) continue;
    sb = si - 7;
    sf = sb + 13;
    ec = (double)(3*c);
    for (tstem = 4; tstem <= 5; tstem++) {
      if (sb >= (sl-8)) goto NX;
      sc = sf;
      sx = si - 2;
      for (tloop = 5; tloop <= 9; tloop++) {
        if (tloop > 7)
          penalty = 3.0*(double)(tloop - tstem - 2);
        else
          penalty = 3.0*(double)(12 - tloop - tstem);
        s1 = sb;
        s2 = sc;
        se = s1 + tstem;
        energy = ec + bem[*se][se[4]] + bem[*s1++][*--s2] - penalty;
        while (s1  < se) energy += bem[*s1++][*--s2];
        energy += G[*sx] + A[sx[1]] + T[sx[3]] + C[sx[4]] + C[sx[5]];
        if (energy >= thresh2)
         { if (i >= nh)
            { fprintf(stderr,"Too many tstem hits\n");
              goto FN; }
           hit[i].pos = sb;
           hit[i].loop = tloop;
           hit[i].stem = tstem;
           hit[i].energy = energy;
           i++; }
        sx++;
        sc++; }
      NX:
      if (--sb < ss) break;
      sf++;
    }
  }

  FN:
  return(i); }

int find_astem5(int *si, int *sl, int *astem3, int n3, trna_loop hit[], int nh, csw *sw){
  int i,k;
  int *s1,*s2,*se;
  unsigned int r,tascanthresh;
  double tastemthresh,energy;
  static unsigned int tem[6] = { 0,0,0,0,0,0 };
  static unsigned int A[6] = { 0,0,0,2,0,0 };
  static unsigned int C[6] = { 0,0,2,0,0,0 };
  static unsigned int G[6] = { 0,2,0,1,0,0 };
  static unsigned int T[6] = { 2,0,1,0,0,0 };
  static double abem[6][6] =
  //    A        C       G       T
   { { -2.144,  -0.428, -2.144, ATBOND, 0.000, 0.000 },
     { -0.428,  -2.144, GCBOND, -2.144, 0.000, 0.000 },
     { -2.144,  GCBOND, -2.144,  1.286, 0.000, 0.000 },
     {  ATBOND, -2.144,  1.286, -0.428, 0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000, 0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000, 0.000, 0.000 } };
  tascanthresh = (unsigned int)sw->tascanthresh;
  tastemthresh = sw->tastemthresh;
  i = 0;
  sl += n3;
  se = astem3 + n3 - 1;
  tem[0] = A[*se];
  tem[1] = C[*se];
  tem[2] = G[*se];
  tem[3] = T[*se];
  while (--se >= astem3)
   { tem[0] = (tem[0] << 4) + A[*se];
     tem[1] = (tem[1] << 4) + C[*se];
     tem[2] = (tem[2] << 4) + G[*se];
     tem[3] = (tem[3] << 4) + T[*se]; }
  r = tem[*si++];
  k = 1;
  while (++k < n3) r = (r >> 4) + tem[*si++];
  while (si < sl)
   { r = (r >> 4) + tem[*si++];
     if ((r & 15) >= tascanthresh)
      { s1 = astem3;
        s2 = si;
        se = s1 + n3;
        energy = abem[*s1++][*--s2];
        while (s1  < se)
         energy += abem[*s1++][*--s2];
        if (energy >= tastemthresh)
         { if (i >= nh)
            { fprintf(stderr,"Too many astem5 hits\n");
              goto FN; }
           hit[i].pos = si - n3;
           hit[i].energy = energy;
           i++; }}}
  FN:
  return(i); }

int aatail(int *s, int *ext, csw *sw){
  int score,e;
  static int A[6] = { 1,0,0,0,0,0 };
  static int C[6] = { 0,1,0,0,0,0 };
  if (sw->aataildiv)
   { score = 0;
     e = 0;
     if (A[s[3]])
      { score++;
        e = 3; }
     if (C[s[2]])
      { score++;
        if (!e) e = 2; }
     if (C[s[1]])
      { score++;
        if (!e) e = 1; }
     if (score < 2)
      if (A[*s]) score++;
     *ext = ++e;
     return(score); }
  else
   { score = 1;
     e = 1;
     if (C[s[1]])
      { score++;
        e = 2;
        if (C[s[2]])
         { score++;
           e = 3;
           if (A[s[3]])
            { score++;
              e = 4; }}}
     *ext = e;
     return(score); }}

void ti_genedetected(int *seq, gene *te, csw *sw) {
  int as,aext,as8,aext8,nbasefext,*s;
  int pseq[2*MAXETRNALEN+1];
  te->nbase = te->astem1 + te->spacer1 + te->spacer2 + 2*te->dstem +
              te->dloop +  2*te->cstem + te->cloop +
              te->var + 2*te->tstem + te->tloop + te->astem2;
  s = te->ps + te->nbase + te->nintron;
  as = aatail(s,&aext,sw);
  if (sw->extastem)
   if (te->astem1 == 7)
    if (bp[te->ps[-1]][*s])
     { as8 = aatail(s+1,&aext8,sw);
       if (as8 >= as)
        { te->ps--;
          te->nbase += 2;
          te->anticodon++;
          if (te->nintron > 0) te->intron++;
          te->astem1 = 8;
          te->astem2 = 8;
          as = as8;
          aext = aext8; }}
  nbasefext = te->nbase + ASTEM2_EXT;
  te->nbase += aext;

  te->start = (long)(te->ps - seq);
  te->stop = (long)(te->ps - seq + te->nbase); }

genes* predict_trnas(dna_sequence *seq, csw *sw) {
  int i,j,k,intron,nt,nth,nd1,nd2,ndx,ndh,na,nah,nppah,nc,nch,tfold,tarm;
  int dstem,dloop,mindist,maxdist,tmindist,tmaxdist,tmmindist,tmmaxdist;
  int tarmthresh,tmstrict,sp2min,sp2max,ige[7];
  int *se,*sc,*sb,*si,*tpos,*tend,*apos,*dpos,*tloopfold,*tmv,*cend;
  int *s1,*s2,*sd,*sf,*sl,*sg1,*sg2,*cposmin,*cposmax,*cpos;
  unsigned int r,q,c;
  double e,ec,he,the,thet,ethresh,energy,cenergy,denergy,ienergy;
  double tdarmthresh,genergy,energy2,energyf,energyf6;

  // linked list of observed genes (e.g., tRNAs)
  genes* gs = NULL;

  static unsigned int TT[6] =
   { 0x00, 0x00, 0x00, 0x11, 0x00, 0x00 };
  static unsigned int GG[6] =
   { 0x00, 0x00, 0x11, 0x00, 0x00, 0x00 };
  static unsigned int ct[6] = { 0,0,0,0,0,0 };
  static unsigned int cA[6] = { 0,0,0,2,0,0 };
  static unsigned int cC[6] = { 0,0,2,0,0,0 };
  static unsigned int cG[6] = { 0,2,0,1,0,0 };
  static unsigned int cT[6] = { 2,0,1,0,0,0 };
  static int yic[9] = { 1,0,0,0,0,0,0,0,0 };
  static int tic[9] = { 1,1,0,0,0,0,0,0,0 };
  static int a1ic[9] = { 1,1,1,0,0,0,0,0,0 };
  static int a2ic[9] = { 1,1,1,1,0,0,0,0,0 };
  static int a3ic[9] = { 1,1,1,1,1,0,0,0,0 };
  static int ric[9] = { 1,1,1,1,1,1,0,0,0 };
  static int goffb[13] = { 0,0,0,0,1,2,2,2,2,2,2,2,2 };
  static int goffe[13] = { 0,0,0,0,2,3,4,4,5,6,6,6,6 };
  static int cY[6] = { 0,1,0,1,0,0 };
  static int cR[6] = { 1,0,1,0,0,0 };
  static double ilw = 0.002;
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,7.0,0.0,0.0 };
  static double Y[6] = { 0.0,3.0,0.0,3.0,0.0,0.0 };
  static double R[6] = { 2.0,0.0,2.0,0.0,0.0,0.0 };
  static double YP[6] = { 0.0,3.0,0.0,3.0,0.0,0.0 };
  static double RP[6] = { 2.0,0.0,2.0,0.0,0.0,0.0 };
  static double RI[6] = { 0.1,0.0,0.05,0.0,0.0,0.0 };
  static double GI[6] = { 0.0,0.0,0.1,0.0,0.0,0.0 };
  static double YI[6] = { 0.0,0.1,0.0,0.1,0.0,0.0 };
  static double AI[6] = { 1.0,0.0,0.0,0.0,0.0,0.0 };
  static double GC[6] = { 0.0,1.5,6.0,0.0,0.0,0.0 };
  static double G3[6] = { 0.0,6.0,12.0,12.0,0.0,0.0 };
  static double dR[6] = { 6.0,0.0,6.0,0.0,0.0,0.0 };
  static double RH[6] = { 3.0,0.0,3.0,0.0,0.0,0.0 };
  static double AGT[6] = { 6.0,0.0,6.0,6.0,0.0,0.0 };
  static double dT[6] = { 0.0,0.0,0.0,6.0,0.0,0.0 };
  static double dbem[6][6] =
  //    A        C       G       T
   { { -2.144,  -0.428, -2.144,  ATBOND, 0.000, 0.000 },
     { -0.428,  -2.144, GCBOND, -2.144,  0.000, 0.000 },
     { -2.144,  GCBOND, -2.144,  1.286,  0.000, 0.000 },
     {  ATBOND, -2.144,  1.286, -0.428,  0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000,  0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000,  0.000, 0.000 } };
  static double dfem[6][6] =
  //    A        C       G       T
   { { -4.000,  -4.000, -4.000,  ATBOND, 0.000, 0.000 },
     { -4.000,  -4.000,  GCBOND, -4.000, 0.000, 0.000 },
     { -4.000,  GCBOND, -4.000,   1.286, 0.000, 0.000 },
     {  ATBOND, -4.000,  1.286,  -4.000, 0.000, 0.000 },
     {  0.000,  0.000,   0.000,   0.000, 0.000, 0.000 },
     {  0.000,  0.000,   0.000,   0.000, 0.000, 0.000 } };
  static double cbem[6][6] =
  //    A            C            G            T
   { { -1.072,      -0.214,      -1.072,       2.0*ATBOND, 0.000, 0.000 },
     { -0.214,      -1.072,       2.0*GCBOND, -1.072,      0.000, 0.000 },
     { -1.072,       2.0*GCBOND, -1.072,       3.400,      0.000, 0.000 },
     {  2.0*ATBOND, -1.072,       3.400,      -0.214,      0.000, 0.000 },
     {  0.000,       0.000,       0.000,       0.000,      0.000, 0.000 },
     {  0.000,       0.000,       0.000,       0.000,      0.000, 0.000 } };

  trna_loop thit[NTH],chit[NC],ahit[NA];
  trna_dloop dhit[ND];

  gene te;

  gene t;
  t.ps        = NULL;
  t.nbase     = 0;
  t.comp      = 0;
  t.start     = 0L;
  t.stop      = 0L;
  t.astem1    = 7;
  t.astem2    = 7;
  t.aatail    = 1;
  t.spacer1   = 2;
  t.spacer2   = 1;
  t.dstem     = 3;
  t.dloop     = 9;
  t.cstem     = 5;
  t.cloop     = 7;
  t.intron    = 0;
  t.nintron   = 0;
  t.anticodon = 0;
  t.var       = 15;
  t.varbp     = 0;
  t.tstem     = 5;
  t.tloop     = 7;
  t.genetype  = tRNA;
  t.energy    = 0.0;
  t.asst      = 0;
  t.tps       = 0;
  t.tpe       = 0;

  ethresh = sw->trnathresh;
  tmmindist = MINTPTSDIST + MINTPDIST;
  tmmaxdist = MAXTPTSDIST + MAXTPDIST;
  tmindist = (MINTRNALEN + sw->minintronlen - MAXTSTEM_DIST);
  tmaxdist = (MAXTRNALEN + sw->maxintronlen - MINTSTEM_DIST);
  mindist = tmindist;
  maxdist = tmaxdist;
  tarmthresh = sw->ttarmthresh;
  tdarmthresh = sw->tdarmthresh;
  tmstrict = sw->tmstrict;
  sp2min = sw->sp2min;
  sp2max = sw->sp2max;

  // effect: updates thit
  nth = find_tstems(seq,thit,NTH,sw);

  nt = -1;
  while (++nt < nth) {
    tpos = thit[nt].pos;
    t.tloop = thit[nt].loop;
    t.tstem = thit[nt].stem;
    tfold = tpos[-1];
    tloopfold = tpos + t.tstem + 1;
    tarm = 2*t.tstem + t.tloop;
    tend = tpos + tarm;
    tmv = tpos - VARMIN;
    te.energy = ethresh;
    the = thit[nt].energy;

    nah = find_astem5(tpos-maxdist,tpos-mindist,tend,7,ahit,NA,sw);
    if (sw->threshlevel < 1.0) { /* find_tstems is generating extra tstems */
       the -= (G[tpos[t.tstem]] + G[tpos[t.tstem+1]]);
       if (the < tarmthresh) continue;
    }
    na = -1;
    while (++na < nah) {
      apos = ahit[na].pos;
      if (apos < (tpos - tmaxdist)) continue;
      if (apos > (tpos - tmindist)) break;
      he = the + ahit[na].energy;

      /* find dstems */
     
      ndh = 0;
      sc = apos + 8;
      energyf = dfem[sc[5]][tfold];
      sl = sc + sw->sp1max;
      while (sc < sl) {
        energy2 = dT[sc[-2]] + RH[*(sc-1)] + GC[*sc] + dfem[sc[-2]][sc[4]];
        energyf6 = dfem[sc[6]][tfold];
        for (dstem = 3; dstem <= 4; dstem++) {
          sd = sc + dstem;
          dloop = 3;
          se = sd + dloop;
          energy = energy2 + 6.0 + dR[*(se-1)] + energyf;
          if (dstem == 3)
           if (energyf < 0.0) energyf = energyf6;
          se += dstem;
          s1 = sc;
          s2 = se;
          sf = s1 + dstem;
          while (s1 < sf) energy += dbem[*s1++][*--s2];
          if (energy >= tdarmthresh) {
            if (ndh >= ND) goto DFL;
            dhit[ndh].pos = sc;
            dhit[ndh].end = se;
            dhit[ndh].loop = dloop;
            dhit[ndh].stem = dstem;
            dhit[ndh].energy = energy;
            ndh++;
          }
          sg1 = sd + 1;
          sg2 = sd + 6;
          q = GG[*sg1++];
          ige[1] = q & 3;
          j = 2;
          while (sg1 <= sg2) {
             q = (q >> 4) + GG[*sg1++];
             ige[j++] = q & 3;
          }
          for (dloop = 4; dloop <= 11; dloop++) {
            j = goffb[dloop];
            k = goffe[dloop];
            c = ige[j++];
            while (j <= k) c = c | ige[j++];
            genergy = G3[c];
            se = sd + dloop;
            energy = energy2 + genergy + dR[*(se-1)] + energyf;
            se += dstem;
            s1 = sc;
            s2 = se;
            sf = s1 + dstem;
            while (s1 < sf) energy += dbem[*s1++][*--s2];
            if (energy >= tdarmthresh) {
               if (ndh >= ND) goto DFL;
               dhit[ndh].pos = sc;
               dhit[ndh].end = se;
               dhit[ndh].loop = dloop;
               dhit[ndh].stem = dstem;
               dhit[ndh].energy = energy;
               ndh++;
            }
          }
        }
        s1 = sc;
        s2 = sc + 16;
        sd = sc + 6;
        j = bp[*s1][*--s2];
        while (++s1 < sd) j += bp[*s1][*--s2];
        if (j >= 6) {
          energy = dT[sc[-1]] + RH[*sc] + GC[*(sc+1)] + energyf6;
          energy += G[*++sd];
          energy += G[*++sd];
          energy += AGT[*++sd] + dfem[sc[-1]][sc[4]];
          sd += 7;
          s1 = sc;
          s2 = sd;
          sf = s1 + 6;
          while (s1 < sf) energy += dbem[*s1++][*--s2];
          if (energy >= tdarmthresh)
           { if (ndh >= ND) goto DFL;
             dhit[ndh].pos = sc;
             dhit[ndh].end = sd;
             dhit[ndh].loop = 4;
             dhit[ndh].stem = 6;
             dhit[ndh].energy = energy;
             ndh++; }
        }
        s1 = sc;
        s2 = sc + 18;
        sd = sc + 7;
        j = bp[*s1][*--s2];
        while (++s1 < sd) j += bp[*s1][*--s2];
        if (j >= 7) {
          energy = energy2 + dfem[sc[7]][tfold];
          energy += G[*++sd];
          energy += G[*++sd];
          energy += AGT[*++sd];
          sd += 8;
          s1 = sc;
          s2 = sd;
          sf = s1 + 7;
          while (s1 < sf) energy += dbem[*s1++][*--s2];
          if (energy >= tdarmthresh)
           { if (ndh >= ND) goto DFL;
             dhit[ndh].pos = sc;
             dhit[ndh].end = sd;
             dhit[ndh].loop = 4;
             dhit[ndh].stem = 7;
             dhit[ndh].energy = energy;
             ndh++; }
        }
        energyf = energyf6;
        sc++;
      }
      goto DFN;
      DFL:
      fprintf(stderr,"Too many D-stem hits\n");
      DFN:
    
      /* End of find dstems routine */
    
      nd1 = ndh;
      while (--nd1 >= 0) {
         dstem = dhit[nd1].stem;
         dpos = dhit[nd1].pos;
         if ((int)(dpos - apos) < 9)
           dhit[nd1].energy -= 3.0;
         if (*tloopfold == Guanine)
          { sb = dpos + dstem + 2;
            sc = sb;
            se = sb + dhit[nd1].loop - 3;
            r = TT[*sb++];
            while (sb < se)
             { r = (r >> 4) + TT[*sb++];
               if (r & 2)
                { dhit[nd1].energy += 10.0;
                  break; }}
            r = GG[*sc++];
            while (sc < se)
             { r = (r >> 4) + GG[*sc++];
               if (r & 2)
               { dhit[nd1].energy -= 12.0;
                 break; }}}
      }
      nd1 = ndh;
      while (--nd1 >= 0) {
        if (!dhit[nd1].end) continue;
        cpos = dhit[nd1].end;
        denergy = dhit[nd1].energy;
        ndx = nd1;
        nd2 = nd1;
        while (--nd2 >= 0) {
          if (dhit[nd2].end != cpos) continue;
          e = dhit[nd2].energy;
          if (e > denergy) {
            denergy = e;
            dhit[ndx].end = NULL;
            ndx = nd2;
          }
        }
      }
      cposmin = 0;
      cposmax = 0;
      nd1 = ndh;
      while (--nd1 >= 0) {
         if (!dhit[nd1].end) continue;
         cposmin = dhit[nd1].end;
         cposmax = cposmin;
         break;
      }
      nd2 = nd1;
      while (--nd2 >= 0) {
        if (!(cpos = dhit[nd2].end)) continue;
        if (cpos < cposmin) cposmin = cpos;
        if (cpos > cposmax) cposmax = cpos;
      }
      for (cpos = cposmin + sp2min; cpos <= (cposmax + sp2max); cpos++) {
        denergy = -INACTIVE;
        ndx = -1;
        nd1 = ndh;
        while (--nd1 >= 0) {
          if (!dhit[nd1].end) continue;
          if ((dhit[nd1].end + sp2max) < cpos) continue;
          if ((dhit[nd1].end + sp2min) > cpos) continue;
          e = dhit[nd1].energy;
          if (e > denergy)
           { denergy = e;
             ndx = nd1; }
        }
        if (ndx < 0) continue;
        denergy += he;
        if (denergy < (te.energy - 49.0)) continue;
    
        /* find cstems */
       
        nch = 0;
        si = cpos;
        sc = cpos + 5;
        se = cpos + 4;
        ct[0] = cA[*se];
        ct[1] = cC[*se];
        ct[2] = cG[*se];
        ct[3] = cT[*se];
     
        while (--se >= cpos) {
          ct[0] = (ct[0] << 4) + cA[*se];
          ct[1] = (ct[1] << 4) + cC[*se];
          ct[2] = (ct[2] << 4) + cG[*se];
          ct[3] = (ct[3] << 4) + cT[*se];
        }
        si += 11;
        se = tmv - VARDIFF - 5;
        if (si < se) si = se;
        r = ct[*si++];
        r = (r >> 4) + ct[*si++];
        r = (r >> 4) + ct[*si++];
        r = (r >> 4) + ct[*si++];
        while (si < tmv) {
          r = (r >> 4) + ct[*si++];
          if ((r & 0xf) >= 5) {
            if (nch >= NC) { fprintf(stderr,"Too many cstem hits\n");
               goto FN;
            }
            chit[nch].pos = si;
            chit[nch].stem = 5;
            chit[nch].loop = (int)(si - sc - 5);
            if (chit[nch].loop == 9)
             if (bp[*sc][si[-6]])
              if (cY[sc[2]])
               if (cR[sc[6]])
                if (cY[sc[1]])
                 { chit[nch].stem = 6;
                   chit[nch].loop = 7; }
            s1 = cpos;
            s2 = si;
            se = s1 + chit[nch].stem;
            chit[nch].energy = cbem[*s1++][*--s2];
            while (s1  < se)
             chit[nch].energy += cbem[*s1++][*--s2];
            nch++; 
          }
        }
        FN:
        
        /* end of find cstems routine */
        
        nc = -1;
        while (++nc < nch) {
          energy = denergy + chit[nc].energy;
          if (energy < (te.energy - 19.0)) continue;
          cend = chit[nc].pos;
          t.var = (int)(tpos - cend);
          t.cloop = chit[nc].loop;
          t.cstem = chit[nc].stem;
          intron = 0;
          if (t.cloop < 9) {
            if (sw->minintronlen > 0) continue;
            if (sw->cloop7)
             if (t.cloop != 7) continue;
            t.nintron = 0;
            if (t.var > 17) energy += vloop_stability(cend,t.var,&t.varbp);
            sb = cpos + t.cstem;
            energy += T[*(sb + 1)] + Y[*(sb)] + R[*(sb + 5)] - 0.05*t.var - ((t.cloop == 7)?0.0:6.0);
          } else {
            t.nintron = t.cloop - 7;
            if (t.nintron > sw->maxintronlen) continue;
            if (t.nintron < sw->minintronlen) continue;
            if (t.var > 17) energy += vloop_stability(cend,t.var,&t.varbp);
            if (energy < (te.energy - 9.0)) continue;
            t.cloop = 7;
            sb = cpos + t.cstem;
            se = sb + t.nintron;
            if (sw->ifixedpos) {
              intron = 6;
              cenergy = YP[*sb] + T[sb[1]] + RP[sb[5]];
            } else {
              cenergy = YP[*se] + T[*(se+1)] + RP[*(se+5)];
              ienergy = cenergy + RI[*sb] + GI[*(se-1)] + AI[se[-2]]*YI[se[-1]];
              for (j = 1; j <= 7; j++) {
                si = se + j - 1;
                ec = YP[*(sb + yic[j]*t.nintron)] + T[*(sb + tic[j]*t.nintron + 1)] + RP[*(sb + ric[j]*t.nintron + 5)];
                e = ec + RI[*(sb + j)] + GI[*si] + AI[si[-1]]*YI[*si];
                if (j == 6) e += 0.01;
                if (e > ienergy) {
                  ienergy = e;
                  cenergy = ec;
                  intron = j;
                }
              }
            }
            energy +=  cenergy - 10.0 - ilw*(t.nintron  + 1.1*t.var);
            if (t.nintron >= 130) {
              si = se + intron;
              j = si[-1];
              if (j != Guanine) {
                 if (si[-2] != Adenine) energy -= 4.0;
                   if (j != Cytosine)
                     if (j != Thymine)
                       energy -= 8.0;
              }
            }
          }

          dstem = dhit[ndx].stem;
          dpos = dhit[ndx].pos;
          if (dstem >= 6) {
            if (sb[2 + a1ic[intron]*t.nintron] != Thymine) continue;
            if (sb[3 + a2ic[intron]*t.nintron] != Cytosine) continue;
            if (sb[4 + a3ic[intron]*t.nintron] != Adenine) continue;
            energy += 3.0;
          } else {
            if (!(dpos[-1] & 5)) {
               i = 0;
               si = cend;
               se = cend + 4;
               while (si < se)
                { if (!(*si++ & 5))
                   { if (++i >= 2)
                      { energy += 3.0;
                        break; }}
                  else
                   i = 0; }
            }
          }
          if (t.cstem >= 6) {
            if (sb[2 + a1ic[intron]*t.nintron] == Cytosine)
             if (sb[3 + a2ic[intron]*t.nintron] == Thymine)
              if (sb[4 + a3ic[intron]*t.nintron] == Adenine)
               energy += 4.0;
          }
          if (energy < ethresh) continue;
          t.energy = energy;
          t.dstem = dstem;
          t.astem1 = (t.dstem < 6)?7:((t.tstem < 5)?9:8);
          t.astem2 = t.astem1;
          t.ps = apos + 7 - t.astem1;
          t.nbase = (int)(tend - t.ps) + t.astem2;
          t.dloop = dhit[ndx].loop;
          t.spacer1 = (int)(dpos - apos) - 7;
          t.spacer2 = (int)(cpos - dhit[ndx].end);
          j = (int)(cpos - t.ps) + t.cstem;
          t.anticodon = j + 2;

          if (t.nintron > 0) {
            t.intron = j + intron;
            if ((t.nbase + t.nintron) > MAXTRNALEN) {
              ti_genedetected(seq->seq,&t,sw);
              gs = add_gene(gs, t);
              continue;
            }
          }
          if (energy < te.energy) continue;
          te = t;
          ti_genedetected(seq->seq,&t,sw);
          gs = add_gene(gs, t);
        }
      }
    }
  }
  return(gs);
}


/* Helpers */

int length(char *s) {
  int i = 0;
  while (*s++) i++;
  return(i); }

dna_sequence* as_sequence(char* dna){
  dna_sequence* d = (dna_sequence*)malloc(sizeof(dna_sequence));
  d->size = length(dna);
  int* seq = (int*)malloc(d->size * sizeof(int));
  for(int i = 0; i < d->size; i++){
    switch(dna[i]) {
      case 'A': 
        seq[i] = Adenine;
        break;
      case 'T': 
        seq[i] = Thymine;
        break;
      case 'G': 
        seq[i] = Guanine;
        break;
      case 'C': 
        seq[i] = Cytosine;
        break;
      default: 
        seq[i] = AMBIG;
    }
  }
  d->seq = seq;
  return(d);  
}


int main(int z, char *v[]) {
  // Default parameters
  static csw sw;

  sw.trna         = 1;
  sw.tmrna        = 0;
  sw.mtrna        = 0;
  sw.cloop7       = 0;
  sw.extastem     = 1;
  sw.aataildiv    = 0;
  sw.sp1max       = 3;
  sw.sp2min       = 0;
  sw.sp2max       = 2;
  sw.maxintronlen = 0;
  sw.minintronlen = 0;
  sw.ifixedpos    = 0;
  sw.tmstrict     = 1;
  sw.loffset      = MAXTAGDIST + 20;
  sw.roffset      = MAXTAGDIST + 20;
  sw.threshlevel  = 1.0;
  sw.trnathresh   = 132.0;
  sw.ttscanthresh = 4.0;
  sw.ttarmthresh  = 29.0;
  sw.tdarmthresh  = 26.0;
  sw.tastemthresh = 7.5;
  sw.tascanthresh = 8.0;

  char* test_seq = "GCAAGATTGCCATCTACTGGGAAGGCGAATTTGTGTGGCAGAACACGGAAGACCAGCGCGACGACATCAAGGGCATCCATCTTCAAGATGATGGCAACTTTGTCTTGTAGTAAGTATCTTAACCTGTGGAGTGGTATTGCAGAAGCTGATTAGTTACTTAGTACCCATGACGACGAAGCTGTTTGGGCTTCTGATACTTGTGGCCAGGGTGACGAGGTCTATCTGGTTGTGCAGGATGATGGAAATGTAGTCTTGTACAAGGGCGAGGATGATGAAGCTGAGGCCGTTTGGGCGACTAGCACCAACCAGATTTAAGCATCGCAAGCTAGTTTTTTTATCTGTGAATAAGAAGGAACATGAAAAGAACAACAGCATTACTACGGCACTCTTATGATTTAAGCCTCAATATGTCCTTTGAATGGTCCAAGCCGTGAATGTCACTGTAAGTAAGTTCAGACTGCATAGTCCCTATTCTGAAGATTGATCAGACATCATCAAGAGCCACGCCCGGTTAGCTCAATCGGTAGAGCGTGAGACTCTTACGAGTACGCGATCTCAAGGTTGCGGGTTCGACCCCCGCATCGGGCTGTTCCTATATTCGATTGCGAGCGTTAGCGAGGCTTTCATTTTTGCATTTTTGTTCTTTGGTACCTTTAAACCTTGCATCTATCTATCTTACGTGCACTGTTCAGTGGATGGAAGAATAAAGTTGTGGACACGTAATGCTCTACTTACTTTGAGTTTATCTTTTTACTTGGAACAGATTTTGAGTTTATGCTCCGCGGACCCCGCCTTGATTTTGTGAGGCTGTGATGTCGCCATTTCCGGAGTGGCATTATCATTTCATCTTTACGCAGATCAACATAAACCACAAAGTTAACCATGCGCATCGAAAAGTATCCCTTATCTCCGCCTAGCTTGGAGGAAGTGGCTGTTAAACTTCAAGCACCCTTGGCAGCCAATTATGAGCACTCCACAGTGAGTGTCGTTTCCTGCCCGGATCTTCGGCAAGCACCCTATCATCTAGCGACCGAGGGGCTATCAGGAGATGAGAAG";

  dna_sequence* seq = as_sequence(test_seq);

  genes* gs = predict_trnas(seq, &sw);

  while(gs && gs->gene){
    fprintf(stdout, "(%d,%d)\n", gs->gene->start, gs->gene->stop);
    gs = gs->next;
  }

  return(0);
}
