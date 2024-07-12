#ifndef __ARAGORN_TRNA_HPP__
#define __ARAGORN_TRNA_HPP__

#include <vector>
#include <string>

typedef struct { long start;
                 long stop;
                 double energy; } hit;

std::vector<hit> predict_trnas(std::string &dna);

#endif
