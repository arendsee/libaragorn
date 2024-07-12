#include "aragorn-trna.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

struct FastaEntry {
    std::string header;
    std::string sequence;
};

std::vector<FastaEntry> read_fasta(const std::string& filename) {
    std::vector<FastaEntry> entries;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return entries;
    }

    std::string line;
    FastaEntry current_entry;
    bool in_sequence = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_sequence) {
                entries.push_back(std::move(current_entry));
                current_entry = FastaEntry();
            }
            in_sequence = true;

            current_entry.header = line.substr(1, line.size() - 1);
        } else {
            current_entry.sequence.append(line);
        }
    }

    if (in_sequence) {
        entries.push_back(std::move(current_entry));
    }

    return entries;
}


int main(int z, char *v[]) {
  std::string filename(v[1]);
  std::vector<FastaEntry> entries = read_fasta(filename);
  for(size_t i = 0; i < entries.size(); i++){
    std::vector<hit> hits = predict_trnas(entries[i].sequence);
    for(size_t j = 0; j < hits.size(); j++){
      fprintf(stdout, "%ld\t%ld\t%ld\t%f\n", i, hits[j].start, hits[j].stop, hits[j].energy);
    }
  }
  return(0);
}
