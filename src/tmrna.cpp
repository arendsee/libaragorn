#include "aragorn.hpp"
#include "common.hpp"

#include <stdio.h>
#include <stdlib.h>

#define INACTIVE        2.0e+35
#define TERM            -1

#define ASTEM2_EXT      9
#define MINTSTEM_DIST   (17 + ASTEM2_EXT)
#define MAXTSTEM_DIST   (26 + ASTEM2_EXT)
#define MAXINTRONLEN    3000
#define MINCTRNALEN     62
#define MAXCTRNALEN     110
#define MINTRNALEN      (MINCTRNALEN + 1)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM2_EXT)
#define MAXETRNALEN     (MAXTRNALEN + MAXINTRONLEN)
#define VARMIN          3
#define MINTPTSDIST     50
#define MAXTPTSDIST     321
#define TPWINDOW        (MAXTPTSDIST - MINTPTSDIST + 1)
#define MINTPDIST       50
#define MAXTPDIST       250
#define TPDISTWINDOW    (MAXTPDIST - MINTPDIST + 1)
#define MINTAGDIST      12
#define MAXTAGDIST      102
#define MINRNACDIST     (MINTPDIST - 5)
#define MAXRNACDIST     (MAXTPDIST - 5)
#define MAXPPINTRONDIST 250
#define TMPTRAILER      145
#define MINPPASDIST     MINTSTEM_DIST
#define MAXPPASDIST     MAXTSTEM_DIST + MAXPPINTRONDIST
#define MINPPTSTPDIST   MINTSTEM_DIST + MINTPDIST
#define MAXPPTSTPDIST   MAXTSTEM_DIST+ASTEM2_EXT+MAXTPDIST+MAXPPINTRONDIST

#define NA          MAXINTRONLEN
#define NH          2000
#define NTH         3000
#define NC          5000
#define ATBOND      2.5

int par_tmstrict = 0;
double par_tmrthresh = 9.0;
int par_aataildiv = 0;
int par_maxintronlen = 0;
int par_minintronlen = 0;
double par_ttarmthresh = 29.0;
double par_tmcathresh = 25.0;
double par_tmrnathresh = 325.0;
double par_tmathresh = 14.0;
double par_tmcthresh = 10.0;
double par_ttscanthresh = 4.0;
double par_tastemthresh = 7.5;
double par_tascanthresh = 8.0;
static int par_loffset = 0;
static int par_roffset = 0;

typedef struct { int *pos;
                 int stem;
                 int loop;
                 double energy; } trna_loop;

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

double find_taghairpin(int *seq)
{ int i,*s,*sb,*se,*sf;
  unsigned int c,m,mx;
  static unsigned int A[6] = { 0,0,0,1,0,0 };
  static unsigned int C[6] = { 0,0,1,0,0,0 };
  static unsigned int G[6] = { 0,1,0,1,0,0 };
  static unsigned int T[6] = { 1,0,1,0,0,0 };
  static unsigned int t[6] = { 0,0,0,0,0,0 };
  mx = 0;
  sb = seq - 20;
  se = seq - 13;
  sf = seq - 4;
  t[0] = A[*sb];
  t[1] = C[*sb];
  t[2] = G[*sb];
  t[3] = T[*sb];
  while (++sb < se)
   { t[0] = (t[0] << 4) | A[*sb];
     t[1] = (t[1] << 4) | C[*sb];
     t[2] = (t[2] << 4) | G[*sb];
     t[3] = (t[3] << 4) | T[*sb]; }
  while (sb < sf)
   { t[0] = ((t[0] << 4) | A[*sb]) & 0xffffffff;
     t[1] = ((t[1] << 4) | C[*sb]) & 0xffffffff;
     t[2] = ((t[2] << 4) | G[*sb]) & 0xffffffff;
     t[3] = ((t[3] << 4) | T[*sb]) & 0xffffffff;
     sb++;
     s = seq + 20;
     se = seq + 2;
     m = t[*s--];
     while (s > se)
      { m = (m >> 4) + t[*s--];
        c = m & 0xf;
        if (c > mx) mx = c; }
     i = 7 - (int)mx;
     while (i-- > 0)
      { m = m >> 4;
        c = m & 0xf;
        if (c > mx) mx = c; }}
  return((double)(mx << 1)); }

double find_tag_upstream_hairpin(int *se)
{ int *sb,*sd,*sf,*sh,*s;
  unsigned int c,m,mx;
  static unsigned int A[6] = { 0,0,0,0x10000,0,0 };
  static unsigned int C[6] = { 0,0,0x10000,0,0,0 };
  static unsigned int G[6] = { 0,0x10000,0,0x10000,0,0 };
  static unsigned int T[6] = { 0x10000,0,0x10000,0,0,0 };
  static unsigned int t[6] = { 0,0,0,0,0,0 };
  mx = 0;
  sf = se - 4;
  sb = se - 20;
  t[0] = A[*se];
  t[1] = C[*se];
  t[2] = G[*se];
  t[3] = T[*se];
  while (--se > sf)
   { t[0] = (t[0] >> 4) | A[*se];
     t[1] = (t[1] >> 4) | C[*se];
     t[2] = (t[2] >> 4) | G[*se];
     t[3] = (t[3] >> 4) | T[*se]; }
  sh = se - 4;
  sd = se - 30;
  while (se > sb)
   { t[0] = ((t[0] >> 4) | A[*se]);
     t[1] = ((t[1] >> 4) | C[*se]);
     t[2] = ((t[2] >> 4) | G[*se]);
     t[3] = ((t[3] >> 4) | T[*se]);
     s = sh;
     m = t[*s];
     while (--s > sd)
       {  m = (m >> 4) + t[*s];
          c = m & 0xf;
          if (c > mx) mx = c;
          if (mx == 5) goto FND; }
     sd--;
     sh--;
     se--; }
  return(0.0);
  FND:
  return(15.0); }

void remove_intron(int *s1, int *s2, int nbase, int intron, int nintron)
{ int *s1e;
  s1e = s1 + intron;
  nbase -= intron;
  while (s1 < s1e) *s2++ = *s1++;
  s1 += nintron;
  s1e = s1 + nbase;
  while (s1 < s1e) *s2++ = *s1++;
  *s2 = TERM; }

int find_resume_seq(int *s, int ls, trna_loop hit[], int nh)
{ int e,i,j,k,a,aa[3],*si,*sb,*sf,*st,*sl;
  double al;
  unsigned int r,c,thresh;
  static int nps[105] =
   { 0,0,0,0, 0,0,0,0,
     0,0,0,0, 1,1,1,1,
     0,0,0,0, 1,1,1,1,
     0,0,0,0, 1,1,1,1,
     0,0,0,0, 1,1,1,1,
     1,1,1,1, 1,1,1,1,
     0,1,0,1, 0,0,0,0,
     0,1,1,1, 1,1,1,1,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,
     0,0,0,0, 0,0,0,0,0 };
  static double score[4] = { 36.0, 66.0, 62.0, 72.0 };
  static unsigned int tem[6] =
   { 0x10310000, 0x01000101, 0x00010030,
     0x02000100, 0x00000000, 0x00000000 };
  static int A[6] = { 0,1,1,1,1,1 };
  static int V[6] = { 0,0,0,1,1,1 };
  static int M[6] = { 0,0,1,1,1,1 };
  thresh = (unsigned int)par_tmrthresh;
  i = 0;
  sl = s + ls;
  r = tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  r = (r >> 4) + tem[*s++];
  if (par_tmstrict)
    while (s < sl)
     { r = (r >> 4) + tem[*s++];
       if ((c = (r & 0xF)) < thresh) continue;
       c -= (V[s[1]] + V[s[2]] + M[s[5]] + A[s[8]]);
       if (c < thresh) continue;
       if (i >= nh) goto FL;
       st = s - 2;
       si = st;
       sb = st + MINTAGDIST + 2;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (*si++ != Thymine)
           si++;
          else
           if (*si == Adenine)
            { if (!(*++si & 5)) goto ST1; }
           else
            if (*si == Guanine)
             { if (*++si == Adenine) goto ST1; }
            else si++;
          si++; }
       continue;
       ST1:
       if (si < sb) continue;
       al = 0.0;
       k = 0;
       j = -11;
       while (j < -2)
     { a = si[j++];
          a = (a << 2) | si[j++];
       if (a == 9) al = (double)(11 + 2*((j + 9)/3));
          a = (a << 2) | si[j++];
          aa[k++] = a; }
       hit[i].pos = st;
       hit[i].stem = (int)(si - st);
       e = (nps[aa[1]] << 1) | (nps[aa[2]]);
       hit[i].energy = (double)(c << 2) + score[e] + al +
                       find_taghairpin(si) +
                       find_tag_upstream_hairpin(st-10);
       i++; }
  else
    while (s < sl)
     { r = (r >> 4) + tem[*s++];
       if ((c = (r & 0xF)) < thresh) continue;
       if (i >= nh) goto FL;
       st = s - 2;
       si = st + MINTAGDIST;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (*si++ != Thymine)
           si++;
          else
           if (*si == Adenine)
            { if (!(*++si & 5)) goto ST2; }
           else
            if (*si == Guanine)
             { if (*++si == Adenine) goto ST2; }
            else si++;
          si++; }
       continue;
       ST2:
       hit[i].pos = st;
       hit[i].stem = (int)(si - st);
       e = (nps[(si[-8] << 4) | (si[-7] << 2) | si[-6]] << 1) |
           (nps[(si[-5] << 4) | (si[-4] << 2) | si[-3]]);
       hit[i].energy = 46.0 + (double)(c << 2) + score[e];
       i++; }
  FN:
  return(i);
  FL:
  fprintf(stderr,"Too many resume sequence hits\n");
  goto FN; }

int *base_copy3(int *from, int *to, int n)
{ while (n-- > 0) *to++ = *from++;
  *to = TERM;
  return(to);  }

tRNA make_tmrna(gene &g){
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

int find_astem5(int *si, int *sl, int *astem3, int n3,
                trna_loop hit[], int nh)
{ int i,k;
  int *s1,*s2,*se;
  unsigned int r,tascanthresh;
  double tastemthresh,energy;
  static unsigned int tem[6] = { 0,0,0,0,0,0 };
  static unsigned int A[6] = { 0,0,0,2,0,0 };
  static unsigned int C[6] = { 0,0,2,0,0,0 };
  static unsigned int G[6] = { 0,2,0,1,0,0 };
  static unsigned int T[6] = { 2,0,1,0,0,0 };
  static double abem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  tascanthresh = (unsigned int)par_tascanthresh;
  tastemthresh = par_tastemthresh;
  i = 0;
  sl += n3;
  se = astem3 + n3 - 1;
  tem[0] = A[*se];
  tem[1] = C[*se];
  tem[2] = G[*se];
  tem[3] = T[*se];
  while (--se >= astem3) {
    tem[0] = (tem[0] << 4) + A[*se];
    tem[1] = (tem[1] << 4) + C[*se];
    tem[2] = (tem[2] << 4) + G[*se];
    tem[3] = (tem[3] << 4) + T[*se];
  }
  r = tem[*si++];
  k = 1;
  while (++k < n3) r = (r >> 4) + tem[*si++];
  while (si < sl) {
    r = (r >> 4) + tem[*si++];
    if ((r & 15) >= tascanthresh) {
      s1 = astem3;
      s2 = si;
      se = s1 + n3;
      energy = abem[*s1++][*--s2];
      while (s1  < se)
       energy += abem[*s1++][*--s2];
      if (energy >= tastemthresh) {
        if (i >= nh) {
          fprintf(stderr,"Too many astem5 hits\n");
          goto FN;
        }
        hit[i].pos = si - n3;
        hit[i].energy = energy;
        i++;
      }
    }
  }
  FN:
  return(i);
}


int find_tstems(int *s, int ls, trna_loop hit[], int nh) {
  int i,r,c,tstem,tloop,ithresh1;
  int *s1,*s2,*se,*ss,*si,*sb,*sc,*sf,*sl,*sx;
  double ec,energy,penalty,thresh2;
  static double bem[6][6] =
   { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
     { -0.428,-2.144, 3.000,-2.144, 0.000, 0.000 },
     { -2.144, 3.000,-2.144, 1.286, 0.000, 0.000 },
     {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static double A[6] = { 2.0,0.0,0.0,0.0,0.0,0.0 };
  static double C[6] = { 0.0,2.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,2.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,2.0,0.0,0.0 };
  static int tem[6] = { 0x0100, 0x0002, 0x2220, 0x0220, 0x0000, 0x0000 };

  i = 0;
  ithresh1 = (int)par_ttscanthresh;
  thresh2 = par_ttarmthresh;
  ss = s + par_loffset;
  si = ss + 4 - 1;
  sl = s + ls - par_roffset + 5 + 3;
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
        if (energy >= thresh2) {
          if (i >= nh) {
            fprintf(stderr,"Too many tstem hits\n");
            goto FN;
          }
          hit[i].pos = sb;
          hit[i].loop = tloop;
          hit[i].stem = tstem;
          hit[i].energy = energy;
          i++;
        }
        sx++;
        sc++;
      }
      NX:
      if (--sb < ss) break;
      sf++;
    }
  }
  FN:
  return(i);
}


void tmopt(std::vector<tRNA> &gs,
          trna_loop *th, int tarm, double the,
          trna_loop *ahit, int nah,
          int *seq) {
  int r,na,nr,nrh,ibase,flag,nbasefext;
  int *s,*v,*s1,*s2,*sa,*sb,*se,*sf,*ps,*tpos,pseq[MAXETRNALEN+1];
  static int gtem[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
  static double A[6] = { 6.0,0.0,0.0,0.0,0.0,0.0 };
  static double Ar[6] = { 10.0,0.0,0.0,0.0,0.0,0.0 };
  static double Cr[6] = { 0.0,10.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double Ga[6] = { 0.0,0.0,7.0,0.0,0.0,0.0 };
  static double K[6] = { 0.0,0.0,6.0,6.0,0.0,0.0 };
  static double Tr[6] = { 0.0,0.0,0.0,10.0,0.0,0.0 };
  double e,energy,penergy,tenergy,aenergy,athresh,cthresh,cathresh;
  static double bem[6][6] =
   { { -1.072,-0.214,-1.072, ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 3.000,-1.072, 0.000, 0.000 },
     { -1.072, 3.000,-1.072, 1.286, 0.000, 0.000 },
     {  ATBOND,-1.072, 1.286,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static trna_loop rhit[NH];
  gene te;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,0,0,0,13,8,0,28,0,0,3,0,5,7,
     1,0.0,0,0,0 };
  tpos = th->pos;
  flag = 0;
  te.energy = par_tmrnathresh;
  athresh = par_tmathresh;
  cthresh = par_tmcthresh;
  cathresh = par_tmcathresh;
  s = tpos + tarm + 4;
  v = tpos + th->stem - 10;
  energy = K[*v] + G[v[1]] + A[v[2]];
  e = K[v[1]] + G[v[2]] + A[v[3]];
  if (e > energy) energy = e;
  if (energy < 18.0) energy = 0.0;
  tenergy = Tr[*s]+Cr[s[1]]+Cr[s[2]]+Ar[s[3]] + energy + 1.59*the;
  nrh = find_resume_seq(tpos-MAXTPTSDIST,TPWINDOW,rhit,NH);
  nr = -1;
  while (++nr < nrh) {
    ps = rhit[nr].pos;
    penergy = tenergy + rhit[nr].energy - 0.001*((double)(tpos - ps));
    if (rhit[nr].stem < 24) penergy -= 15.0;
    na = -1;
    while (++na < nah) {
      aenergy = ahit[na].energy;
      if (aenergy < athresh) continue;
      t.ps = ahit[na].pos;
      if (t.ps < (ps - MAXTPDIST)) continue;
      if (t.ps > (ps - MINTPDIST)) break;
      energy = -INACTIVE;
      sa = t.ps + t.astem1;
      for (sb=sa+9, se=sb+t.cstem; sb <= (sa+16); sb++,se++)
      for (sf = tpos-3; sf >= (tpos-7); sf--) {
        s1 = sb;
        s2 = sf;
        e = bem[*s1++][*--s2];
        while (s1 < se) e += bem[*s1++][*--s2];
        if (e > energy) {
          energy = e;
          t.var = (int)(tpos - sf);
          t.dloop = (int)(sb - sa);
        }
      }
      if (energy < cthresh) continue;
      energy += aenergy;
      if (energy < cathresh) continue;
      sb = sa + 3;
      sf = sa + 7;
      r = gtem[*sb++];
      while (sb < sf) {
        r = (r >> 4) + gtem[*sb++];
        if ((r & 3) == 2) {
          energy += 14.0;
          break;
        }
      }
      t.energy = penergy + Ga[t.ps[1]] + Ga[t.ps[2]] + energy;
      if (t.energy > te.energy) {
        flag = 1;
        t.tstem = th->stem;
        t.tloop = th->loop;
        t.tps = (int)(ps - t.ps);
        t.tpe = t.tps + rhit[nr].stem;
        ibase = (int)(tpos - t.ps);
        t.nintron = ibase - t.var - 2*t.cstem - t.dloop - t.astem1;
        t.nbase = ibase + tarm + t.astem2 - t.nintron;
        te = t;
      }
    }
  }
  if (flag){
     te.start = (long)(te.ps - seq);
     s = te.ps + te.nbase + te.nintron;
     nbasefext = te.nbase + ASTEM2_EXT;

     te.intron = te.astem1 + te.dloop + te.cstem;
     te.asst = 0;
     base_copy3(te.ps,te.eseq,nbasefext+te.nintron);
     remove_intron(te.ps,pseq,nbasefext,
                   te.intron,te.nintron);
     base_copy3(pseq,te.seq,te.nbase);
     gs.push_back(make_tmrna(te));
  }
}

void tmopt_perm(std::vector<tRNA> &gs,
          trna_loop *th, int tarm, double the,
          trna_loop *ahit, int nah,
          int *seq)
{ int r,na,nr,nrh,flag;
  int *s,*v,*s1,*s2,*sa,*sb,*se,*sf,*ps,*apos,*tpos;
  static int gtem[6] = { 0x00,0x00,0x11,0x00,0x00,0x00 };
  double e,energy,penergy,tenergy,aenergy,athresh,cthresh,cathresh;
  static double A[6] = { 6.0,0.0,0.0,0.0,0.0,0.0 };
  static double Ar[6] = { 10.0,0.0,0.0,0.0,0.0,0.0 };
  static double Cr[6] = { 0.0,10.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static double Ga[6] = { 0.0,0.0,7.0,0.0,0.0,0.0 };
  static double K[6] = { 0.0,0.0,6.0,6.0,0.0,0.0 };
  static double Tr[6] = { 0.0,0.0,0.0,10.0,0.0,0.0 };
  static double bem[6][6] =
   { { -1.072,-0.214,-1.072, ATBOND, 0.000, 0.000 },
     { -0.214,-1.072, 3.000,-1.072, 0.000, 0.000 },
     { -1.072, 3.000,-1.072, 1.286, 0.000, 0.000 },
     {  ATBOND,-1.072, 1.286,-0.214, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
     {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };
  static trna_loop rhit[NH];
  gene te;
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,0,0,0,13,8,0,28,0,0,3,0,5,7,
     1,0.0,0,0,0 };
  tpos = th->pos;
  flag = 0;
  te.energy = par_tmrnathresh;
  athresh = par_tmathresh;
  cthresh = par_tmcthresh;
  cathresh = par_tmcathresh;
  s = tpos + tarm + 4;
  v = tpos + th->stem - 10;
  energy = K[*v] + G[v[1]] + A[v[2]];
  e = K[v[1]] + G[v[2]] + A[v[3]];
  if (e > energy) energy = e;
  if (energy < 18.0) energy = 0.0;
  tenergy = Tr[*s]+Cr[s[1]]+Cr[s[2]]+Ar[s[3]]+ energy + 1.59*the;
  na = -1;
  while (++na < nah){
    aenergy = ahit[na].energy;
    if (aenergy < athresh) continue;
    apos = ahit[na].pos;
    if (apos < (tpos + MINTSTEM_DIST)) continue;
    if (apos > (tpos + MAXTSTEM_DIST + MAXPPINTRONDIST)) break;
    energy = -INACTIVE;
    sa = apos + t.astem1;
    for (sb=sa+9, se=sb+t.cstem; sb <= (sa+16); sb++,se++) {
      for (sf = tpos-3; sf >= (tpos-7); sf--) {
        s1 = sb;
        s2 = sf;
        e = bem[*s1++][*--s2];
        while (s1 < se) e += bem[*s1++][*--s2];
        if (e > energy) {
          energy = e;
          t.var = (int)(tpos - sf);
          t.dloop = (int)(sb - sa);
        }
      }
    }
    if (energy < cthresh) continue;
    energy += aenergy;
    if (energy < cathresh) continue;
    sb = sa + 3;
    sf = sa + 7;
    r = gtem[*sb++];
    while (sb < sf) {
      r = (r >> 4) + gtem[*sb++];
      if ((r & 3) == 2) {
        energy += 14.0;
        break;
      }
    }
    penergy = tenergy + Ga[apos[1]] + Ga[apos[2]] + energy;
    nrh = find_resume_seq(apos+MINTPDIST,TPWINDOW,rhit,NH);
    nr = -1;
    while (++nr < nrh) {
      ps = rhit[nr].pos;
      t.energy = penergy + rhit[nr].energy;
      if (rhit[nr].stem < 24) t.energy -= 15.0;
      if (t.energy > te.energy) {
        flag = 1;
        t.tstem = th->stem;
        t.tloop = th->loop;
        t.asst = (long)(apos - tpos) + t.var + t.cstem;
        t.ps = tpos - t.var - t.cstem;
        t.tps = (int)(ps - t.ps);
        t.tpe = t.tps + rhit[nr].stem;
        te = t;
      }
    }
  }
  if (flag) {
    te.start = (long)(te.ps - seq) - 54;
    te.intron = te.cstem + te.var + 2*te.tstem + te.tloop +
                te.astem2;
    base_copy3(te.ps-54,te.eseq,te.tpe+1+TMPTRAILER);
    te.nbase = te.astem1 + te.dloop + te.cstem;
    base_copy3(te.ps+te.asst,te.seq,te.nbase);
    base_copy3(te.ps,te.seq+te.nbase,te.intron + ASTEM2_EXT);
    te.nbase += te.intron;
    te.nintron = te.tpe - te.nbase + 1 + TMPTRAILER;
    te.intron += 54;
    te.tps += 54;
    te.tpe += 54;
    te.asst += 54;

    gs.push_back(make_tmrna(te));
  }
}


std::vector<tRNA> predict_tmrnas(std::string &dna) {

  std::vector<int> seq_vec = dna2int(dna);
  int* seq = &seq_vec[0];
  int lseq = seq_vec.size();

  std::vector<tRNA> gs;

  int nt,nth,nah,nppah,tarm;
  int mindist,maxdist,tmindist,tmaxdist,tmmindist,tmmaxdist;
  int *tpos,*tend;
  double the,thet;
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  static trna_loop thit[NTH],ahit[NA];
  static gene t =
   { "",{TERM},{TERM},NULL,0,0,0L,0L,7,7,1,2,1,3,9,5,7,0,0,0,15,0,5,7,
     1,0.0,0,0,0 };
  tmmindist = MINTPTSDIST + MINTPDIST;
  tmmaxdist = MAXTPTSDIST + MAXTPDIST;
  tmindist = (MINTRNALEN + par_minintronlen - MAXTSTEM_DIST);
  tmaxdist = (MAXTRNALEN + par_maxintronlen - MINTSTEM_DIST);
  mindist = (tmindist < tmmindist)?tmindist:tmmindist;
  maxdist = (tmaxdist > tmmaxdist)?tmaxdist:tmmaxdist;
  nth = find_tstems(seq,lseq,thit,NTH);
  nt = -1;
  while (++nt < nth) {
    tpos = thit[nt].pos;
    t.tloop = thit[nt].loop;
    t.tstem = thit[nt].stem;
    tarm = 2*t.tstem + t.tloop;
    tend = tpos + tarm;
    the = thit[nt].energy;
    nah = find_astem5(tpos-maxdist,tpos-mindist,tend,7,ahit,NA);
    thet = the - G[tpos[t.tstem]] - G[tpos[t.tstem+1]];
    if (par_tmstrict){
      if (thet >= par_ttarmthresh)
        tmopt(gs, thit+nt,tarm,thet,ahit,nah,seq);
    } else {
      tmopt(gs, thit+nt,tarm,the,ahit,nah,seq);
    }
    nppah = find_astem5(tpos+MINPPASDIST,tpos+MAXPPASDIST, tend,7,ahit+nah,NA-nah);
    tmopt_perm(gs, thit+nt,tarm,the,ahit+nah,nppah,seq);
    if (thet < par_ttarmthresh) continue;
    the = thet;
  }
  return(gs);
}
