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

double find_taghairpin(const std::vector<int> &seq, int seq_idx)
{ int i, s, sb, se, sf;
  unsigned int c,m,mx;
  static unsigned int A[6] = { 0,0,0,1,0,0 };
  static unsigned int C[6] = { 0,0,1,0,0,0 };
  static unsigned int G[6] = { 0,1,0,1,0,0 };
  static unsigned int T[6] = { 1,0,1,0,0,0 };
  static unsigned int t[6] = { 0,0,0,0,0,0 };
  mx = 0;
  sb = seq_idx - 20;
  se = seq_idx - 13;
  sf = seq_idx - 4;
  t[0] = A[seq[sb]];
  t[1] = C[seq[sb]];
  t[2] = G[seq[sb]];
  t[3] = T[seq[sb]];
  while (++sb < se)
   { t[0] = (t[0] << 4) | A[seq[sb]];
     t[1] = (t[1] << 4) | C[seq[sb]];
     t[2] = (t[2] << 4) | G[seq[sb]];
     t[3] = (t[3] << 4) | T[seq[sb]]; }
  while (sb < sf)
   { t[0] = ((t[0] << 4) | A[seq[sb]]) & 0xffffffff;
     t[1] = ((t[1] << 4) | C[seq[sb]]) & 0xffffffff;
     t[2] = ((t[2] << 4) | G[seq[sb]]) & 0xffffffff;
     t[3] = ((t[3] << 4) | T[seq[sb]]) & 0xffffffff;
     sb++;
     s = seq_idx + 20;
     se = seq_idx + 2;
     m = t[seq[s--]];
     while (s > se)
      { m = (m >> 4) + t[seq[s--]];
        c = m & 0xf;
        if (c > mx) mx = c; }
     i = 7 - (int)mx;
     while (i-- > 0)
      { m = m >> 4;
        c = m & 0xf;
        if (c > mx) mx = c; }}
  return((double)(mx << 1)); }

double find_tag_upstream_hairpin(const std::vector<int> &seq, int se)
{ int sb, sd, sf, sh, s;
  unsigned int c, m, mx;
  static unsigned int A[6] = { 0,0,0,0x10000,0,0 };
  static unsigned int C[6] = { 0,0,0x10000,0,0,0 };
  static unsigned int G[6] = { 0,0x10000,0,0x10000,0,0 };
  static unsigned int T[6] = { 0x10000,0,0x10000,0,0,0 };
  static unsigned int t[6] = { 0,0,0,0,0,0 };
  mx = 0;
  sf = se - 4;
  sb = se - 20;
  t[0] = A[seq[se]];
  t[1] = C[seq[se]];
  t[2] = G[seq[se]];
  t[3] = T[seq[se]];
  while (--se > sf)
   { t[0] = (t[0] >> 4) | A[seq[se]];
     t[1] = (t[1] >> 4) | C[seq[se]];
     t[2] = (t[2] >> 4) | G[seq[se]];
     t[3] = (t[3] >> 4) | T[seq[se]]; }
  sh = se - 4;
  sd = se - 30;
  while (se > sb)
   { t[0] = ((t[0] >> 4) | A[seq[se]]);
     t[1] = ((t[1] >> 4) | C[seq[se]]);
     t[2] = ((t[2] >> 4) | G[seq[se]]);
     t[3] = ((t[3] >> 4) | T[seq[se]]);
     s = sh;
     m = t[seq[s]];
     while (--s > sd)
       {  m = (m >> 4) + t[seq[s]];
          c = m & 0xf;
          if (c > mx) mx = c;
          if (mx == 5) goto FND; }
     sd--;
     sh--;
     se--; }
  return(0.0);
  FND:
  return(15.0); }

std::vector<TrnaLoop> find_resume_seq(const std::vector<int> &seq, int seq_idx, int window_size)
{ int e,i,j,k,a,aa[3],si,sb,sf,st,sl;
  double al;
  unsigned int r,c,thresh;

  std::vector<TrnaLoop> hits;
  int hit_pos;
  int hit_stem;
  double hit_energy;

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
  sl = seq_idx + window_size;
  r = tem[seq[seq_idx++]];
  r = (r >> 4) + tem[seq[seq_idx++]];
  r = (r >> 4) + tem[seq[seq_idx++]];
  r = (r >> 4) + tem[seq[seq_idx++]];
  r = (r >> 4) + tem[seq[seq_idx++]];
  r = (r >> 4) + tem[seq[seq_idx++]];
  r = (r >> 4) + tem[seq[seq_idx++]];
  if (par_tmstrict)
    while (seq_idx < sl)
     { r = (r >> 4) + tem[seq[seq_idx++]];
       if ((c = (r & 0xF)) < thresh) continue;
       c -= (V[seq[seq_idx + 1]] + V[seq[seq_idx + 2]] + M[seq[seq_idx + 5]] + A[seq[seq_idx + 8]]);
       if (c < thresh) continue;
       st = seq_idx - 2;
       si = st;
       sb = st + MINTAGDIST + 2;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (seq[si++] != Thymine)
           si++;
          else
           if (seq[si] == Adenine)
            { if (!(seq[++si] & 5)) goto ST1; }
           else
            if (seq[si] == Guanine)
             { if (seq[++si] == Adenine) goto ST1; }
            else si++;
          si++; }
       continue;
       ST1:
       if (si < sb) continue;
       al = 0.0;
       k = 0;
       j = -11;
       while (j < -2) {
         a = seq[si+j++];
         a = (a << 2) | seq[si+j++];
         if (a == 9) al = (double)(11 + 2*((j + 9)/3));
         a = (a << 2) | seq[si+j++];
         aa[k++] = a;
       }

       hit_pos = st;
       hit_stem = (int)(si - st);
       e = (nps[aa[1]] << 1) | (nps[aa[2]]);
       hit_energy = (double)(c << 2) + score[e] + al +
                    find_taghairpin(seq, si) +
                    find_tag_upstream_hairpin(seq, st-10);

       hits.push_back(make_trna_loop(hit_pos, -1, hit_stem, hit_energy));
     }
  else
    while (seq_idx < sl) {
      r = (r >> 4) + tem[seq[seq_idx++]];
       if ((c = (r & 0xF)) < thresh) continue;
       st = seq_idx - 2;
       si = st + MINTAGDIST;
       sf = st + MAXTAGDIST;
       while (si < sf)
        { if (seq[si++] != Thymine)
           si++;
          else
           if (seq[si] == Adenine)
            { if (!(seq[++si] & 5)) goto ST2; }
           else
            if (seq[si] == Guanine)
             { if (seq[++si] == Adenine) goto ST2; }
            else si++;
          si++; }
       continue;
       ST2:

       hit_pos = st;
       hit_stem = (int)(si - st);
       e = (nps[(seq[si-8] << 4) | (seq[si-7] << 2) | seq[si-6]] << 1) |
           (nps[(seq[si-5] << 4) | (seq[si-4] << 2) | seq[si-3]]);
       hit_energy = 46.0 + (double)(c << 2) + score[e];

       hits.push_back(make_trna_loop(hit_pos, -1, hit_stem, hit_energy));
    }

  return(hits);
}

void tmopt(const std::vector<int> &seq,
           std::vector<tRNA> &gs,
           TrnaLoop tloop, int tarm, double the,
           const std::vector<TrnaAstem> &astem_hits) {
  int r,ibase,flag,nbasefext;
  int s, v, s1, s2, sa, sb, se, sf, ps, tpos, pseq[MAXETRNALEN+1];
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
  Gene te = Gene();
  Gene t = Gene();
  t.spacer1 = 0;
  t.spacer2 = 0;
  t.dstem = 0;
  t.dloop = 13;
  t.cstem = 8;
  t.cloop = 0;
  t.intron = 28;
  t.var = 3;

  tpos = tloop.pos;
  flag = 0;
  te.energy = par_tmrnathresh;
  athresh = par_tmathresh;
  cthresh = par_tmcthresh;
  cathresh = par_tmcathresh;
  s = tpos + tarm + 4;
  v = tpos + tloop.stem - 10;
  energy = K[seq[v]] + G[seq[v+1]] + A[seq[v+2]];
  e = K[seq[v+1]] + G[seq[v+2]] + A[seq[v+3]];
  if (e > energy) energy = e;
  if (energy < 18.0) energy = 0.0;
  tenergy = Tr[seq[s]]+Cr[seq[s+1]]+Cr[seq[s+2]]+Ar[seq[s+3]] + energy + 1.59*the;

  std::vector<TrnaLoop> rhits = find_resume_seq(seq, tpos-MAXTPTSDIST,TPWINDOW);

  for(auto & rhit : rhits) {
    ps = rhit.pos;
    penergy = tenergy + rhit.energy - 0.001*((double)(tpos - ps));
    if (rhit.stem < 24) penergy -= 15.0;

    for(auto & astem : astem_hits) {
      aenergy = astem.energy;
      if (aenergy < athresh) continue;
      t.ps = astem.pos;
      if (t.ps < (ps - MAXTPDIST)) continue;
      if (t.ps > (ps - MINTPDIST)) break;
      energy = -INACTIVE;
      sa = t.ps + t.astem1;
      for (sb=sa+9, se=sb+t.cstem; sb <= (sa+16); sb++,se++)
      for (sf = tpos-3; sf >= (tpos-7); sf--) {
        s1 = sb;
        s2 = sf;
        e = bem[seq[s1++]][seq[--s2]];
        while (s1 < se) e += bem[seq[s1++]][seq[--s2]];
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
      r = gtem[seq[sb++]];
      while (sb < sf) {
        r = (r >> 4) + gtem[seq[sb++]];
        if ((r & 3) == 2) {
          energy += 14.0;
          break;
        }
      }
      t.energy = penergy + Ga[seq[t.ps+1]] + Ga[seq[t.ps+2]] + energy;
      if (t.energy > te.energy) {
        flag = 1;
        t.tstem = tloop.stem;
        t.tloop = tloop.loop;
        t.tps = (int)(ps - t.ps);
        t.tpe = t.tps + rhit.stem;
        ibase = (int)(tpos - t.ps);
        t.nintron = ibase - t.var - 2*t.cstem - t.dloop - t.astem1;
        t.nbase = ibase + tarm + t.astem2 - t.nintron;
        te = t;
      }
    }
  }
  if (flag){
     te.start = (long)(te.ps);
     s = te.ps + te.nbase + te.nintron;
     nbasefext = te.nbase + ASTEM2_EXT;

     te.intron = te.astem1 + te.dloop + te.cstem;
     te.asst = 0;
     gs.push_back(make_trna(te));
  }
}

void tmopt_perm(const std::vector<int> &seq,
          std::vector<tRNA> &gs,
          TrnaLoop tloop, int tarm, double the,
          const std::vector<TrnaAstem> &astem_hits)
{ int r,na,nr,nrh,flag;
  int s, v, s1, s2, sa, sb, se, sf, ps, apos, tpos;
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
  Gene te = Gene();
  Gene t = Gene();
  t.spacer1 = 0;
  t.spacer2 = 0;
  t.dstem = 0;
  t.dloop = 13;
  t.dstem = 0;
  t.dloop = 13;
  t.cstem = 8;
  t.cloop = 0;
  t.intron = 28;
  t.nintron = 0;
  t.var = 3;

  tpos = tloop.pos;
  flag = 0;
  te.energy = par_tmrnathresh;
  athresh = par_tmathresh;
  cthresh = par_tmcthresh;
  cathresh = par_tmcathresh;
  s = tpos + tarm + 4;
  v = tpos + tloop.stem - 10;
  energy = K[seq[v]] + G[seq[v+1]] + A[seq[v+2]];
  e = K[seq[v+1]] + G[seq[v+2]] + A[seq[v+3]];
  if (e > energy) energy = e;
  if (energy < 18.0) energy = 0.0;
  tenergy = Tr[seq[s]]+Cr[seq[s+1]]+Cr[seq[s+2]]+Ar[seq[s+3]]+ energy + 1.59*the;
  for(auto & astem : astem_hits) {
    aenergy = astem.energy;
    if (aenergy < athresh) continue;
    apos = astem.pos;
    if (apos < (tpos + MINTSTEM_DIST)) continue;
    if (apos > (tpos + MAXTSTEM_DIST + MAXPPINTRONDIST)) break;
    energy = -INACTIVE;
    sa = apos + t.astem1;
    for (sb=sa+9, se=sb+t.cstem; sb <= (sa+16); sb++,se++) {
      for (sf = tpos-3; sf >= (tpos-7); sf--) {
        s1 = sb;
        s2 = sf;
        e = bem[seq[s1++]][seq[--s2]];
        while (s1 < se) e += bem[seq[s1++]][seq[--s2]];
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
    r = gtem[seq[sb++]];
    while (sb < sf) {
      r = (r >> 4) + gtem[seq[sb++]];
      if ((r & 3) == 2) {
        energy += 14.0;
        break;
      }
    }
    penergy = tenergy + Ga[seq[apos+1]] + Ga[seq[apos+2]] + energy;
    std::vector<TrnaLoop> rhits = find_resume_seq(seq, apos + MINTPDIST, TPWINDOW);
    for(auto & rhit : rhits) {
      ps = rhit.pos;
      t.energy = penergy + rhit.energy;
      if (rhit.stem < 24) t.energy -= 15.0;
      if (t.energy > te.energy) {
        flag = 1;
        t.tstem = tloop.stem;
        t.tloop = tloop.loop;
        t.asst = apos - tpos + t.var + t.cstem;
        t.ps = tpos - t.var - t.cstem;
        t.tps = (int)(ps - t.ps);
        t.tpe = t.tps + rhit.stem;
        te = t;
      }
    }
  }
  if (flag) {
    te.start = te.ps - 54;
    te.intron = te.cstem + te.var + 2*te.tstem + te.tloop +
                te.astem2;
    te.nbase = te.astem1 + te.dloop + te.cstem;
    te.nbase += te.intron;
    te.nintron = te.tpe - te.nbase + 1 + TMPTRAILER;
    te.intron += 54;
    te.tps += 54;
    te.tpe += 54;
    te.asst += 54;

    gs.push_back(make_trna(te));
  }
}


std::vector<tRNA> predict_tmrnas(std::string &dna) {

  std::vector<int> seq = dna2int(dna);

  std::vector<tRNA> gs;
  std::vector<TrnaAstem> astem_hits, astem_hits_2;

  int nt,nth,nah,nppah,tarm;
  int mindist,maxdist,tmindist,tmaxdist,tmmindist,tmmaxdist;
  int tpos, tend;
  double the,thet;
  static double G[6] = { 0.0,0.0,6.0,0.0,0.0,0.0 };
  Gene t = Gene();

  tmmindist = MINTPTSDIST + MINTPDIST;
  tmmaxdist = MAXTPTSDIST + MAXTPDIST;
  tmindist = (MINTRNALEN + par_minintronlen - MAXTSTEM_DIST);
  tmaxdist = (MAXTRNALEN + par_maxintronlen - MINTSTEM_DIST);
  mindist = (tmindist < tmmindist)?tmindist:tmmindist;
  maxdist = (tmaxdist > tmmaxdist)?tmaxdist:tmmaxdist;

  for(auto & tloop : find_tstems(seq, par_loffset, par_roffset, par_ttscanthresh, par_ttarmthresh)){
    tpos = tloop.pos;
    t.tloop = tloop.loop;
    t.tstem = tloop.stem;
    tarm = 2*t.tstem + t.tloop;
    tend = tpos + tarm;
    the = tloop.energy;
    astem_hits = find_astem5(seq, tpos - maxdist, tpos - mindist, tend, 7, par_tascanthresh, par_tastemthresh);
    thet = the - G[seq[tpos+t.tstem]] - G[seq[tpos+t.tstem+1]];
    if (par_tmstrict){
      if (thet >= par_ttarmthresh)
        tmopt(seq, gs, tloop, tarm, thet, astem_hits);
    } else {
      tmopt(seq, gs, tloop, tarm, the, astem_hits);
    }
    astem_hits_2 = find_astem5(seq, tpos + MINPPASDIST, tpos + MAXPPASDIST, tend, 7, par_tascanthresh, par_tastemthresh);
    tmopt_perm(seq, gs, tloop, tarm, the, astem_hits_2);
    if (thet < par_ttarmthresh) continue;
    the = thet;
  }

  return(gs);
}
