#include "aragorn.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

class FastaEntry {
public:
    std::string header;
    std::string sequence;

    FastaEntry(const std::string& h, const std::string& s) : header(h), sequence(s) {}
};




std::vector<FastaEntry> read_fasta(const std::string& filename) {
    std::vector<FastaEntry> entries;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return entries;
    }

    std::string line;
    FastaEntry current_entry = FastaEntry("", "");
    bool in_sequence = false;

    auto trim = [](std::string& s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    };

    while (std::getline(file, line)) {
        trim(line);
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_sequence && current_entry.sequence.size() > 0) {
                entries.push_back(current_entry);
            } else {
                in_sequence = true;
            }

            current_entry = FastaEntry(line.substr(1), "");
        } else {
            current_entry.sequence += line;
        }
    }

    if (in_sequence && current_entry.sequence.size() > 0) {
        entries.push_back(current_entry);
    }

    return entries;
}



int main(int z, char *v[]) {
  std::string filename(v[1]);
  std::vector<FastaEntry> entries = read_fasta(filename);
  for(size_t i = 0; i < entries.size(); i++){
    std::vector<tRNA> hits = predict_trnas(entries[i].sequence);
    for(size_t j = 0; j < hits.size(); j++){
      // It seems that their start and stop indices are incorrect?
      std::string trna_seq = entries[i].sequence.substr(hits[j].start, hits[j].stop - hits[j].start);
      std::string anticodon = entries[i].sequence.substr(hits[j].start + hits[j].anticodon, 3);
      int component_length = hits[j].astem1 +
                             hits[j].spacer1 +
                             hits[j].dstem * 2 +
                             hits[j].dloop +
                             hits[j].spacer2 +
                             hits[j].cstem * 2 +
                             hits[j].cloop +
                             hits[j].var +
                             hits[j].tstem * 2 +
                             hits[j].tloop +
                             hits[j].astem2;

      std::cout << entries[i].header << std::endl
                << "seq = " << entries[i].sequence << std::endl
                << "tRNA = " << trna_seq << std::endl
                << "anticodon = " << anticodon << std::endl
                << "start = " << hits[j].start << std::endl
                << "stop = " << hits[j].stop << std::endl
                << "length = (" << component_length << "," << hits[j].stop - hits[j].start << ")" << std::endl
                << "intron = " << hits[j].intron << std::endl
                << "nintron = " << hits[j].nintron << std::endl
                << "energy = " << hits[j].energy << std::endl
                << "----------------------------" << std::endl
                << "astem1 = " << hits[j].astem1 << std::endl
                << "spacer1 = " << hits[j].spacer1 << std::endl
                << "dstem = " << hits[j].dstem << std::endl
                << "dloop = " << hits[j].dloop << std::endl
                << "spacer2 = " << hits[j].spacer2 << std::endl
                << "cstem = " << hits[j].cstem << std::endl
                << "cloop = " << hits[j].cloop << std::endl
                << "var = " << hits[j].var << std::endl
                << "tstem = " << hits[j].tstem << std::endl
                << "tloop = " << hits[j].tloop << std::endl
                << "astem2 = " << hits[j].astem2 << std::endl
                << "//" << std::endl;

    }
  }
  return(0);
}

/*

8  8  2  1  4  9   5  7  0  0  35  5   0  5  7  157.108
aaGCGGAGTTAGTTTAGTCTGGTATGACGTCAGCTTCCCAAGCTGAAGGCCGCGGGTTCAAATCCCGCACTCCGCAt

aaGCGGAGTTAGTTTAGTCTGGTATGACGTCAGCTTCCCAAGCTGAAGGCCGCGGGTTCAAATCCCGCACTCCGCAt
  -------                                                           -------
  CGCCTCA

                            cca   3' amino-acyl acceptor
                    5'    a       NCCA is the canonical motif
                       g-c
                       c-g
                       g-c  A-stem
                       g-c
                       a-t
                       g-c
             space1    t-a     ta
                      t   cgccc  a
              tga    a    !!!!!  a   T-loop
             c   tttg     gcggg  c
     D-loop  t   :+!!    c     tt
             g   tgac     c
              gta    g     g
                      t-aag   var
            spacer2   c-g
                      a-t
                      g-c
                      c-g
                     t   a
                     t   a   C-loop (anticodon stem-loop)
                      ccc
                      \ /
                       anticodon


[agt]cca
1797 a
1408 acca
1185 t
867 ac
724 gcca
668 g
345 gc
266 tc
239 tcca
137 c
111 cc
91 acc
54 ccca
34 gcc
20 tcc
18 ccc



seq = catggctcaaGCGGAGTTAGTTTAGTCTGGTATGACGTCAGCTTCCCAAGCTGAAGGCCGCGGGTTCAAATCCCGCACTCCGCAtcaaatttctcta
tRNA = GCGGAGTTAGTTTAGTCTGGTATGACGTCAGCTTCCCAAGCTGAAGGCCGCGGGTTCAAATCCCGCACTCCGCA
start = 9
stop = 85

astem1 = 8
spacer1 = 2
dstem = 4
dloop = 9
spacer2 = 1
cstem = 5
cloop = 7
var = 5
tstem = 5
tloop = 7
astem2 = 8

intron = 0
nintron = 0
anticodon = 35
energy = 157.108
*/
