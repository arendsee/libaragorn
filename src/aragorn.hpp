#ifndef __ARAGORN_HPP__
#define __ARAGORN_HPP__

#include <vector>
#include <string>

typedef struct {
    int start;
    int stop;
    int anticodon;
    double score;
    int astem1;
    int spacer1;
    int dstem;
    int dloop;
    int spacer2;
    int cstem;
    int cloop;
    int intron_start;
    int intron_length;
    int var;
    int tstem;
    int tloop;
    int astem2;
} tRNA;

std::vector<tRNA> predict_trnas(std::string &dna);

std::vector<tRNA> predict_mtrnas(std::string &dna);

std::vector<tRNA> predict_tmrnas(std::string &dna);

#endif
