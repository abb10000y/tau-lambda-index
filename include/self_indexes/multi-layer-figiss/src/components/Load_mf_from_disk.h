#pragma once
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>

class Load_mf_from_disk
{
private:
    std::vector<std::pair<size_t, size_t>> mfSet;
    std::string path;
    void fetch_min_factors();
public:
    Load_mf_from_disk(std::string filePath);    
    std::vector<std::pair<size_t, size_t>> get_min_factors();
};

Load_mf_from_disk::Load_mf_from_disk(std::string filePath) {
    path = filePath;
    fetch_min_factors();
}

void Load_mf_from_disk::fetch_min_factors() {
    std::ifstream ff (path);
    std::string tmp;
    for (int i = 0; i < 15; i++) { getline(ff, tmp); } // depends on the text format
    while (getline(ff, tmp)) {
        for (char& c : tmp) {
            if (!isdigit(c)) { c = ' '; }
        }
        std::stringstream stm(tmp);
        std::size_t a, b;
        stm >> a >> b;
        mfSet.push_back({a, b});
    }
    ff.close();     
}

std::vector<std::pair<size_t, size_t>> Load_mf_from_disk::get_min_factors() { return mfSet; }