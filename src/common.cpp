#include <algorithm>
#include "common.hpp"

std::vector<int> dna2int (const std::string& dna){
  std::vector<int> seq(dna.size());
  for(size_t i = 0; i < dna.size(); i++){
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
  return seq;
}

struct CompareByStart {
    bool operator()(const tRNA& a, const tRNA& b) const {
        return a.start < b.start;
    }
};

std::vector<tRNA> best_hit(std::vector<tRNA> hits) {
  // sort hits by start position
  std::sort(hits.begin(), hits.end(), CompareByStart());

  std::vector<tRNA> best_hits;

  // initialize
  if (hits.size() > 0) {
    best_hits.push_back(hits[0]);
  }

  // for hits that overlap, keep only the highest energy one
  tRNA last;
  for(auto & hit : hits){
    last = best_hits[best_hits.size() - 1];
    if (hit.start > last.stop) {
      best_hits.push_back(hit);
    } else if (hit.score > last.score) {
      best_hits[best_hits.size() - 1] = hit;
    }
  }

  return best_hits;
}

std::vector<TrnaLoop> find_tstems(const std::vector<int>& s,
                                  int loffset, int roffset, double ttscanthresh, double ttarmthresh) {
  std::vector<TrnaLoop> hits;

  int r,c,tstem,tloop;
  int s1, s2, se, si, sb, sc, sf, sl, sx;
  double ec,energy,penalty;
  static double bem[6][6] = {
  //    A        C       G       T
     { -2.144,  -0.428, -2.144,  ATBOND, 0.000, 0.000 }, // A
     { -0.428,  -2.144,  GCBOND, -2.144, 0.000, 0.000 }, // C
     { -2.144,   GCBOND, -2.144,  1.286, 0.000, 0.000 }, // G
     {  ATBOND, -2.144,  1.286, -0.428,  0.000, 0.000 }, // T
     {  0.000,   0.000,  0.000,  0.000,  0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000,  0.000, 0.000 } };

  static double A[6] = { 2.0,0.0,0.0,0.0,0.0,0.0 };
  static double C[6] = { 0.0,2.0,0.0,0.0,0.0,0.0 };
  static double G[6] = { 0.0,0.0,2.0,0.0,0.0,0.0 };
  static double T[6] = { 0.0,0.0,0.0,2.0,0.0,0.0 };
                      // A       C       G       T       -       -
  static int tem[6]  = { 0x0100, 0x0002, 0x2000, 0x0220, 0x0000, 0x0000 };

  // pointer to the moving t-start?
  // And 4 - 1 because???
  si = loffset + 4 - 1;

  // right ending position after applying offset
  // And 5 + 3 because???
  sl = s.size() - roffset + 5 + 3;

  // initialize the scanning bit pattern with the first three bases
  r = tem[s[si++]];
  r = (r >> 4) + tem[s[si++]];
  r = (r >> 4) + tem[s[si++]];

  while (si < sl) {
    r = (r >> 4) + tem[s[si++]];

    // The `r & 0xF` trick gets the rightmost byte from r
    // The 4th byte contains information about the 4 prior bytes
    // that match the pattern GTTC. The match can contain any two of these
    // bases (at the default ttscanthresh and with tem_trna).
    if ((c = (r & 0xF)) < (int)ttscanthresh) continue;

    // T-Loop with canonical numbering:
    //
    //                           (60) (59)
    //  (65) (64) (63) (62) (61)          (58)
    //    |    |    |    |    |            (57)
    //  (49) (50) (51) (52) (53)          (56)
    //                       G   (54) (55) C
    //                            T    T   *
    //
    // subtract 7 from the start index to move from position 56 at the end of
    // the GTTC match to the start of the loop
    sb = si - 7;

    // add 13 to the end index (the extra bases allow for variation in loop size
    sf = sb + 13;

    ec = (double)(3*c);

    // Loop through possible lengths of the T-stem
    for (tstem = 4; tstem <= 5; tstem++) {
      if (sb < (sl-8)){
        sc = sf;
        sx = si - 2;
        // Loop through possible sizes of the T-loop
        for (tloop = 5; tloop <= 9; tloop++) {
          if (tloop > 7)
            penalty = 3.0*(double)(tloop - tstem - 2);
          else
            penalty = 3.0*(double)(12 - tloop - tstem);
          s1 = sb;
          s2 = sc;
          se = s1 + tstem;
          energy = ec + bem[s[se]][s[se + 4]] + bem[s[s1++]][s[--s2]] - penalty;
          while (s1 < se) energy += bem[s[s1++]][s[--s2]];
          energy += G[s[sx]] + A[s[sx + 1]] + T[s[sx + 3]] + C[s[sx + 4]] + C[s[sx + 5]];
          if (energy >= ttarmthresh) {
             hits.push_back(make_trna_loop(sb, tloop, tstem, energy));
          }
          sx++;
          sc++;
        }
      }
      if (--sb < loffset) break;
      sf++;
    }
  }
  return(hits);
}

std::vector<TrnaAstem> find_astem5(const std::vector<int>& seq, int si, int sl, int astem3, int n3, double tascanthresh, double tastemthresh){

  // fprintf(stderr, "entering find_astem5 (si=%d, sl=%d, astem3=%d, n3=%d, seq.size()=%d)\n", si, sl, astem3, n3, seq.size());

  std::vector<TrnaAstem> hits;
  int k;
  int s1, s2, se;
  int r;
  int N = (int) seq.size();
  double energy;
  static int tem[6] = { 0,0,0,0,0,0 };
  static int A[6] = { 0,0,0,2,0,0 };
  static int C[6] = { 0,0,2,0,0,0 };
  static int G[6] = { 0,2,0,1,0,0 };
  static int T[6] = { 2,0,1,0,0,0 };
  static double abem[6][6] =
  //    A        C       G       T
   { { -2.144,  -0.428, -2.144, ATBOND, 0.000, 0.000 },
     { -0.428,  -2.144, GCBOND, -2.144, 0.000, 0.000 },
     { -2.144,  GCBOND, -2.144,  1.286, 0.000, 0.000 },
     {  ATBOND, -2.144,  1.286, -0.428, 0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000, 0.000, 0.000 },
     {  0.000,   0.000,  0.000,  0.000, 0.000, 0.000 } };
  sl += n3;
  se = astem3 + n3 - 1;

  // fprintf(stderr, " > sl=%d se=%d\n", sl, se);

  if(se >= N || si >= N) return hits;
  tem[0] = A[seq[se]];
  tem[1] = C[seq[se]];
  tem[2] = G[seq[se]];
  tem[3] = T[seq[se]];
  while (--se >= astem3)
   { tem[0] = (tem[0] << 4) + A[seq[se]];
     tem[1] = (tem[1] << 4) + C[seq[se]];
     tem[2] = (tem[2] << 4) + G[seq[se]];
     tem[3] = (tem[3] << 4) + T[seq[se]]; }
  r = tem[seq[si++]];
  k = 1;

  // fprintf(stderr, " > si=%d se=%d\n", si, se);

  while (++k < n3) {
    if(si >= N) return hits;
    r = (r >> 4) + tem[seq[si++]];
  }

  // fprintf(stderr, " > si=%d\n", si);

  if(sl >= N) return hits;
  while (si < sl){
    r = (r >> 4) + tem[seq[si++]];
    if ((r & 15) >= (int)tascanthresh){
      s1 = astem3;
      s2 = si;
      se = s1 + n3;
      energy = abem[seq[s1++]][seq[--s2]];
      while (s1 < se) {
        energy += abem[seq[s1++]][seq[--s2]];
      }
      if (energy >= tastemthresh) {
        hits.push_back(make_astem(si - n3, energy));
      }
    }
  }

  // fprintf(stderr, " > si=%d s1=%d s2=%d\n", si, s1, s2);

  return(hits);
}

TrnaLoop make_trna_loop(int pos, int loop, int stem, double energy){
  TrnaLoop x;
  x.pos = pos;
  x.loop = loop;
  x.stem = stem;
  x.energy = energy;
  return(x);
}

TrnaAstem make_astem(int pos, double energy){
  TrnaAstem x;
  x.pos = pos;
  x.energy = energy;
  return(x);
}

tRNA make_trna(Gene &g){
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
    h.tps = g.tps;
    h.tpe = g.tpe;
    return h;
}
