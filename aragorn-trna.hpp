#ifndef __ARAGORN_TRNA_HPP__
#define __ARAGORN_TRNA_HPP__

#include <vector>
#include <string>

typedef struct {
    long start;
    long stop;
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
    int varbp;
    int tstem;
    int tloop;
    double energy;
} tRNA;

std::vector<tRNA> predict_trnas(std::string &dna);

#endif
