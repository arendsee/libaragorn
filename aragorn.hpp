#ifndef __ARAGORN_TRNA_HPP__
#define __ARAGORN_TRNA_HPP__

#include <vector>
#include <string>

typedef struct {
    int start;
    int stop;
    int anticodon;
    int astem1;
    int spacer1;
    int dstem;
    int dloop;
    int spacer2;
    int cstem;
    int cloop;
    int intron;
    int nintron;
    int var;
    int tstem;
    int tloop;
    int astem2;
    double energy;
} tRNA;

std::vector<tRNA> predict_trnas(std::string &dna);

#endif
