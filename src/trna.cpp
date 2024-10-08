#include "aragorn.hpp"
#include "common.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#define INACTIVE        2.0e+35

#define MINTSTEM_DIST   (17 + ASTEM2_EXT)
#define MAXTSTEM_DIST   (26 + ASTEM2_EXT)
#define MINCTRNALEN     62
#define MAXCTRNALEN     110
#define MINTRNALEN      (MINCTRNALEN + 1)
#define MAXTRNALEN      (MAXCTRNALEN + ASTEM2_EXT)
#define VARMIN          3
#define VARDIFF         23                /* VARMAX - VARMIN */
#define MAXTAGDIST      102

#define NC          5000

double bem[6][6] =
 { { -2.144,-0.428,-2.144, ATBOND, 0.000, 0.000 },
   { -0.428,-2.144, GCBOND,-2.144, 0.000, 0.000 },
   { -2.144, GCBOND,-2.144, 1.286, 0.000, 0.000 },
   {  ATBOND,-2.144, 1.286,-0.428, 0.000, 0.000 },
   {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 },
   {  0.000, 0.000, 0.000, 0.000, 0.000, 0.000 } };


class TrnaDloop {
  public:
    int pos;
    int end;
    int stem;
    int loop;
    double energy;
};

class Config {
  public:
    int cloop7;
    int sp1max;
    int sp2min;
    int sp2max;
    int maxintronlen;
    int minintronlen;
    int ifixedpos;
    int loffset;
    int roffset;
    double threshlevel;
    double trnathresh;
    double ttscanthresh;
    double ttarmthresh;
    double tdarmthresh;
    double tastemthresh;
    double tascanthresh;

    Config()
    {
      cloop7       = 1;
      sp1max       = 3;
      sp2min       = 0;
      sp2max       = 2;
      maxintronlen = 300;
      minintronlen = 0;
      ifixedpos    = 0;
      loffset      = 0;
      roffset      = 20;
      threshlevel  = 1.0;
      trnathresh   = 132.0;
      ttscanthresh = 4.0;
      ttarmthresh  = 29.0;
      tdarmthresh  = 26.0;
      tastemthresh = 7.5;
      tascanthresh = 8.0;
    }

};

/* LIBRARY */

double vloop_stability(const std::vector<int>& seq, int sb, int var){
  int e,stem,vstem,loop;
  int sn, sen, se, sc, sd, sf, s;
  int c,cn,m;
  static int A[6] = { 0,0,0x100,0x400,0,0 };
  static int C[6] = { 0,0,0x400,0,0,0 };
  static int G[6] = { 0x100,0x400,0,0x200,0,0 };
  static int T[6] = { 0x400,0,0x200,0,0,0 };
  static int te[6] = { 0,0,0,0,0,0 };
  e = 0;
  sc = sb + 3;
  se = sb + var - 2;
  sf = se - 2;
  te[0] = A[seq[se]];
  te[1] = C[seq[se]];
  te[2] = G[seq[se]];
  te[3] = T[seq[se]];
  while (--se > sf)
   { te[0] = (te[0] >> 4) | A[seq[se]];
     te[1] = (te[1] >> 4) | C[seq[se]];
     te[2] = (te[2] >> 4) | G[seq[se]];
     te[3] = (te[3] >> 4) | T[seq[se]]; }
  while (se >= sc)
   { te[0] = ((te[0] >> 4) | A[seq[se]]);
     te[1] = ((te[1] >> 4) | C[seq[se]]);
     te[2] = ((te[2] >> 4) | G[seq[se]]);
     te[3] = ((te[3] >> 4) | T[seq[se]]);
     s = se - 5;
     sd = se - 7;
     m = te[seq[s]];
     while (--s > sd) m = (m >> 4) + te[seq[s]];
     while (s >= sb)
       {  m = (m >> 4) + te[seq[s]];
          c = m & 0xf;
          if (c >= 9)
           { stem = 3;
             loop = se - s - 3;
             sen = se;
             sn = s + 2;
             while (loop >= 6)
              { if ((cn = vbp[seq[sen - 1]][seq[sn + 1]]) <= 0) break;
                c += cn;
                stem++;
                loop -= 2;
                sen--;
                sn++; }
             if (c > e)
              { e = c;
                vstem = stem; }}
          s--; }
      se--; }
  if (e > 0) {
    return((double)(3*(vstem - 4)));
  } else {
     return(-12.0);
  }
}

void ti_genedetected(Gene& te, const Config& sw) {
  te.nbase = te.astem1 + te.spacer1 + te.spacer2 + 2*te.dstem +
              te.dloop +  2*te.cstem + te.cloop +
              te.var + 2*te.tstem + te.tloop + te.astem2;
  te.start = te.ps;
  te.stop = te.ps + te.nbase; }

std::vector<tRNA> predict_trnas(std::string &dna) {

  std::vector<int> seq = dna2int(dna);

  int i,j,k,intron,nd1,nd2,ndx,ndh,nc,nch,tfold,tarm;
  int dstem,dloop,tmindist,tmaxdist;
  int ige[7];

  int se, sc, sb, si, tpos, tend, apos, dpos, tloopfold, tmv, cend;
  int s1, s2, sd, sf, sl, sg1, sg2, cposmin, cposmax, cpos;

  int r,q,c;
  double e,ec,he,the,energy,cenergy,denergy,ienergy;
  double genergy,energy2,energyf,energyf6;

  std::vector<tRNA> gs;

  static int TT[6] = { 0x00, 0x00, 0x00, 0x11, 0x00, 0x00 };
  static int GG[6] = { 0x00, 0x00, 0x11, 0x00, 0x00, 0x00 };
  static int ct[6] = { 0,0,0,0,0,0 };
  static int cA[6] = { 0,0,0,2,0,0 };
  static int cC[6] = { 0,0,2,0,0,0 };
  static int cG[6] = { 0,2,0,1,0,0 };
  static int cT[6] = { 2,0,1,0,0,0 };
  static int yic[9]  = { 1,0,0,0,0,0,0,0,0 };
  static int tic[9]  = { 1,1,0,0,0,0,0,0,0 };
  static int a1ic[9] = { 1,1,1,0,0,0,0,0,0 };
  static int a2ic[9] = { 1,1,1,1,0,0,0,0,0 };
  static int a3ic[9] = { 1,1,1,1,1,0,0,0,0 };
  static int ric[9]  = { 1,1,1,1,1,1,0,0,0 };
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

  TrnaLoop chit[NC];
  TrnaDloop dhit[ND];

  static Config sw = Config();

  Gene te = Gene();
  Gene t = Gene();

  int N = (int) seq.size();

  tmindist = (MINTRNALEN + sw.minintronlen - MAXTSTEM_DIST);
  tmaxdist = (MAXTRNALEN + sw.maxintronlen - MINTSTEM_DIST);

  // Try to make a tRNA gene for each predicted T-loop
  for(auto & tloop : find_tstems(seq, sw.loffset, sw.roffset, sw.ttscanthresh, sw.ttarmthresh))
  {
    tpos = tloop.pos;
    t.tloop = tloop.loop;
    t.tstem = tloop.stem;
    tfold = seq[tpos - 1];
    tloopfold = tpos + t.tstem + 1;
    tarm = 2*t.tstem + t.tloop;
    tend = tpos + tarm;

    tmv = tpos - VARMIN;
    te.energy = sw.trnathresh;
    the = tloop.energy; 

    // abort if the observed T-loop is too low confidence
    if (sw.threshlevel < 1.0) {
       the -= (G[seq[tpos + t.tstem]] + G[seq[tpos + t.tstem + 1]]);
       if (the < sw.ttarmthresh) continue;
    }

    int astem_max_start = tpos - tmaxdist;
    astem_max_start = astem_max_start < 0 ? 0 : astem_max_start;

    // Look for A-stems in compatible with the current T-loop
    for(auto & astem5 : find_astem5(seq, astem_max_start, tpos-tmindist, tend, 7, sw.tascanthresh, sw.tastemthresh))
    {

      apos = astem5.pos;

      if (apos < (tpos - tmaxdist)) continue;
      if (apos > (tpos - tmindist)) break;

      he = the + astem5.energy;

      /* find dstems */
     
      ndh = 0;
      sc = apos + 8;
      energyf = dfem[seq[sc+5]][tfold];

      sl = sc + sw.sp1max;
      if (sl >= N) break;

      while (sc < sl) {
        energy2 = dT[seq[sc - 2]] + RH[seq[sc - 1]] + GC[seq[sc]] + dfem[seq[sc - 2]][seq[sc + 4]];
        energyf6 = dfem[seq[sc+6]][tfold];

        for (dstem = 3; dstem <= 4; dstem++) {
          sd = sc + dstem;
          if (sd >= N) break;

          dloop = 3;
          se = sd + dloop;
          if (se >= N) break;

          energy = energy2 + 6.0 + dR[seq[se - 1]] + energyf;
          if (dstem == 3)
           if (energyf < 0.0) energyf = energyf6;
          se += dstem;
          if (se >= N) break;

          s1 = sc;
          s2 = se;
          sf = s1 + dstem;
          if (sf >= N) break;

          while (s1 < sf) energy += dbem[seq[s1++]][seq[--s2]];
          if (energy >= sw.tdarmthresh) {
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
          if (sg2 >= N) break;

          q = GG[seq[sg1++]];
          ige[1] = q & 3;
          j = 2;
          while (sg1 <= sg2) {
             q = (q >> 4) + GG[seq[sg1++]];
             ige[j++] = q & 3;
          }
          for (dloop = 4; dloop <= 11; dloop++) {
            j = goffb[dloop];
            k = goffe[dloop];
            c = ige[j++];
            while (j <= k) c = c | ige[j++];
            genergy = G3[c];
            se = sd + dloop;
            if (se >= N) break;

            energy = energy2 + genergy + dR[seq[se-1]] + energyf;
            se += dstem;
            s1 = sc;
            s2 = se;
            sf = s1 + dstem;
            if (se >= N) break;

            while (s1 < sf) energy += dbem[seq[s1++]][seq[--s2]];
            if (energy >= sw.tdarmthresh) {
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
        if (s2 >= N) break;

        j = bp[seq[s1]][seq[--s2]];
        while (++s1 < sd) j += bp[seq[s1]][seq[--s2]];
        if (j >= 6) {
          energy = dT[seq[sc - 1]] + RH[seq[sc]] + GC[seq[sc + 1]] + energyf6;
          if ((sd + 1) >= N) break;
          energy += G[seq[++sd]];
          energy += G[seq[++sd]];
          energy += AGT[seq[++sd]] + dfem[seq[sc - 1]][seq[sc + 4]];
          sd += 7;
          s1 = sc;
          s2 = sd;
          sf = s1 + 6;
          if (sf >= N) break;
          while (s1 < sf) energy += dbem[seq[s1++]][seq[--s2]];
          if (energy >= sw.tdarmthresh)
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
        j = bp[seq[s1]][seq[--s2]];
        while (++s1 < sd) j += bp[seq[s1]][seq[--s2]];
        if (j >= 7) {
          energy = energy2 + dfem[seq[sc+7]][tfold];
          energy += G[seq[++sd]];
          energy += G[seq[++sd]];
          energy += AGT[seq[++sd]];
          sd += 8;
          s1 = sc;
          s2 = sd;
          sf = s1 + 7;
          while (s1 < sf) energy += dbem[seq[s1++]][seq[--s2]];
          if (energy >= sw.tdarmthresh)
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
         if (seq[tloopfold] == Guanine)
          { sb = dpos + dstem + 2;
            sc = sb;
            se = sb + dhit[nd1].loop - 3;
            r = TT[seq[sb++]];
            while (sb < se) {
              r = (r >> 4) + TT[seq[sb++]];
              if (r & 2) {
                 dhit[nd1].energy += 10.0;
                 break;
              }
            }
            r = GG[seq[sc++]];
            while (sc < se)
             { r = (r >> 4) + GG[seq[sc++]];
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
            dhit[ndx].end = 0;
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
      for (cpos = cposmin + sw.sp2min; cpos <= (cposmax + sw.sp2max); cpos++) {
        denergy = -INACTIVE;
        ndx = -1;
        nd1 = ndh;
        while (--nd1 >= 0) {
          if (!dhit[nd1].end) continue;
          if ((dhit[nd1].end + sw.sp2max) < cpos) continue;
          if ((dhit[nd1].end + sw.sp2min) > cpos) continue;
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
        ct[0] = cA[seq[se]];
        ct[1] = cC[seq[se]];
        ct[2] = cG[seq[se]];
        ct[3] = cT[seq[se]];
     
        while (--se >= cpos) {
          ct[0] = (ct[0] << 4) + cA[seq[se]];
          ct[1] = (ct[1] << 4) + cC[seq[se]];
          ct[2] = (ct[2] << 4) + cG[seq[se]];
          ct[3] = (ct[3] << 4) + cT[seq[se]];
        }
        si += 11;
        se = tmv - VARDIFF - 5;
        if (si < se) si = se;
        r = ct[seq[si++]];
        r = (r >> 4) + ct[seq[si++]];
        r = (r >> 4) + ct[seq[si++]];
        r = (r >> 4) + ct[seq[si++]];
        while (si < tmv) {
          r = (r >> 4) + ct[seq[si++]];
          if ((r & 0xf) >= 5) {
            if (nch >= NC) { fprintf(stderr,"Too many cstem hits\n");
               goto FN;
            }
            chit[nch].pos = si;
            chit[nch].stem = 5;
            chit[nch].loop = (int)(si - sc - 5);
            if (chit[nch].loop == 9)
             if (bp[seq[sc]][seq[si - 6]])
              if (cY[seq[sc+2]])
               if (cR[seq[sc+6]])
                if (cY[seq[sc+1]])
                 { chit[nch].stem = 6;
                   chit[nch].loop = 7; }
            s1 = cpos;
            s2 = si;
            se = s1 + chit[nch].stem;
            chit[nch].energy = cbem[seq[s1++]][seq[--s2]];
            while (s1  < se)
             chit[nch].energy += cbem[seq[s1++]][seq[--s2]];
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
            if (sw.minintronlen > 0) continue;
            if (sw.cloop7 && t.cloop != 7) {
              continue;
            }
            t.nintron = 0;
            if (t.var > 17) energy += vloop_stability(seq, cend, t.var);
            sb = cpos + t.cstem;
            energy += T[seq[sb + 1]] + Y[seq[sb]] + R[seq[sb + 5]] - 0.05*t.var - ((t.cloop == 7)?0.0:6.0);
          } else {
            t.nintron = t.cloop - 7;
            if (t.nintron > sw.maxintronlen) continue;
            if (t.nintron < sw.minintronlen) continue;
            if (t.var > 17) energy += vloop_stability(seq, cend, t.var);
            if (energy < (te.energy - 9.0)) continue;
            t.cloop = 7;
            sb = cpos + t.cstem;
            se = sb + t.nintron;
            if (sw.ifixedpos) {
              intron = 6;
              cenergy = YP[seq[sb]] + T[seq[sb+1]] + RP[seq[sb+5]];
            } else {
              cenergy = YP[seq[se]] + T[seq[se+1]] + RP[seq[se+5]];
              ienergy = cenergy + RI[seq[sb]] + GI[seq[se-1]] + AI[seq[se - 2]]*YI[seq[se - 1]];
              for (j = 1; j <= 7; j++) {
                si = se + j - 1;
                ec = YP[seq[sb + yic[j]*t.nintron]] + T[seq[sb + tic[j]*t.nintron + 1]] + RP[seq[sb + ric[j]*t.nintron + 5]];
                e = ec + RI[seq[sb + j]] + GI[seq[si]] + AI[seq[si - 1]]*YI[seq[si]];
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
              j = seq[si - 1];
              if (j != Guanine) {
                 if (seq[si - 2] != Adenine) energy -= 4.0;
                   if (j != Cytosine)
                     if (j != Thymine)
                       energy -= 8.0;
              }
            }
          }

          dstem = dhit[ndx].stem;
          dpos = dhit[ndx].pos;
          if (dstem >= 6) {
            if (seq[sb + 2 + a1ic[intron] * t.nintron] != Thymine) continue;
            if (seq[sb + 3 + a2ic[intron] * t.nintron] != Cytosine) continue;
            if (seq[sb + 4 + a3ic[intron] * t.nintron] != Adenine) continue;
            energy += 3.0;
          } else {
            if (!(seq[dpos - 1] & 5)) {
               i = 0;
               si = cend;
               se = cend + 4;
               while (si < se)
                { if (!(seq[si++] & 5))
                   { if (++i >= 2)
                      { energy += 3.0;
                        break; }}
                  else
                   i = 0; }
            }
          }
          if (t.cstem >= 6) {
            if (seq[sb + 2 + a1ic[intron]*t.nintron] == Cytosine)
             if (seq[sb + 3 + a2ic[intron]*t.nintron] == Thymine)
              if (seq[sb + 4 + a3ic[intron]*t.nintron] == Adenine)
               energy += 4.0;
          }
          if (energy < sw.trnathresh) continue;
          t.energy = energy;
          t.dstem = dstem;

          // What is this rule for? Why should A-stem length depend on d and t
          // stem lengths? Shouldn't it depend on the length of the base-pairing
          // with A-stem 2?
          t.astem1 = (t.dstem < 6)?7:((t.tstem < 5)?9:8);

          // astem2 length is always just set to the length of astem1
          t.astem2 = t.astem1;
          t.ps = apos + 7 - t.astem1;
          t.nbase = tend - t.ps + t.astem2;
          t.dloop = dhit[ndx].loop;
          t.spacer1 = dpos - apos - 7;
          t.spacer2 = cpos - dhit[ndx].end;
          j = cpos - t.ps + t.cstem;
          t.anticodon = j + 2;

          if (t.nintron > 0) {
            t.intron = j + intron;
            if ((t.nbase + t.nintron) > MAXTRNALEN) {
              ti_genedetected(t,sw);
              gs.push_back(make_trna(t));
              continue;
            }
          }
          if (energy < te.energy) continue;
          te = t;
          ti_genedetected(t,sw);
          gs.push_back(make_trna(t));
        }
      }
    }
  }

  gs = best_hit(gs);

  return(gs);
}
