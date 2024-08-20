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
