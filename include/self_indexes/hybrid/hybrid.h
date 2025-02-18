#pragma once

#include "hsi_parser_lz77/compute_lz77.hpp"
#include "index/HybridSelfIndex.h"

#ifndef uchar
#define uchar unsigned char
#endif

#ifndef ulong
#define ulong unsigned long
#endif

class hybrid {
private:
    HybridSelfIndex *index;
    void LZ_parsing(const std::string &rawText);

    template<typename char_type, typename text_offset_type>
    void store_lz77(std::vector<std::pair<text_offset_type, text_offset_type> > &parsing, ulong nT, char *fileStore) {
    uint64_t n_phrases = parsing.size();
    uint64_t pos, len, MAX_len, MAX_pos, i;
    uint lg_n, lg_len;

    fprintf(stderr, "\nTo Store Dictionary:\n");
    MAX_len = MAX_pos = 0;
    for (i = 0; i < n_phrases; ++i) {
        pos = parsing[i].first;
        len = parsing[i].second;
        //fprintf(stderr, "(%lu,%lu) ", pos, len);
        if (pos > MAX_pos)
            MAX_pos = pos;
        if (len > MAX_len)
            MAX_len = len;
    }
    lg_n = 1 + (uint)(log(MAX_pos)/log(2));
    lg_len = 1 + (uint)(log(MAX_len)/log(2));

    fprintf(stderr, "\nlg_n = %d, MAX_len = %lu, lg_len = %d \n", lg_n, MAX_len, lg_len);

    ulong size_ARR_POS = n_phrases*lg_n/(8*sizeof(ulong));		// number of (ulong) cells for ARR_POS
    if ((n_phrases*lg_n)%(8*sizeof(ulong)))
        size_ARR_POS++;
    ulong *ARR_POS = new ulong[size_ARR_POS];

    ulong size_ARR_LEN = n_phrases*lg_len/(8*sizeof(ulong));
    if ((n_phrases*lg_len)%(8*sizeof(ulong)))
        size_ARR_LEN++;
    ulong *ARR_LEN = new ulong[size_ARR_LEN];

    fprintf(stderr, "\n\nStoring Dictionary of %lu phrases\n", n_phrases);

    for (i=0; i < n_phrases; ++i) {
        pos = parsing[i].first;
        len = parsing[i].second;
        setNum64(ARR_POS, (ulong)i*lg_n, lg_n, pos);
        setNum64(ARR_LEN, (ulong)i*lg_len, lg_len, len);
        //fprintf(stderr, "(%lu, %lu) ", pos, len);
    }

    fprintf(stderr, "\nSave data structure in file %s\n", fileStore);

    fprintf(stderr, "\n n_phrases = %lu\n", n_phrases);
    std::ofstream os(fileStore, std::ofstream::binary);
    os.write((const char*)&nT, sizeof(ulong));
    os.write((const char*)&n_phrases, sizeof(uint64_t));	// number of factors
    os.write((const char*)&lg_n, sizeof(uint));
    os.write((const char*)&lg_len, sizeof(uint));
    os.write((const char*)ARR_POS, size_ARR_POS*sizeof(ulong));
    os.write((const char*)ARR_LEN, size_ARR_LEN*sizeof(ulong));
    os.close();
    }
    
public:
    hybrid(){}
    ~hybrid(){ index->~HybridSelfIndex(); }
    void serialized(){ index->saveStructure(); }
    void build(char* inputFilePath, uint M, char* outputFolderPath);
    void load(const std::string &inputFolderPath);
    void locate(uchar *pat, uint m, ulong *nOcc, ulong **occ);
};

void hybrid::LZ_parsing(const std::string &rawTextPath) {
    typedef uint8_t char_type;
    typedef uint64_t text_offset_type;
    typedef std::pair<text_offset_type, text_offset_type> pair_type;

    std::ifstream file(rawTextPath, std::ios::binary | std::ios::ate);
    if (!file) {
        std::cerr << "Error: Cannot open file " << rawTextPath << std::endl;
        return;
    }

    text_offset_type text_length = file.tellg();
    file.seekg(0, std::ios::beg);

    char_type* text = new char_type[text_length];
    file.read(reinterpret_cast<char*>(text), text_length);
    file.close();
    

    std::vector<pair_type> parsing;
    compute_lz77(text, text_length, parsing);

    // Store to disk
    std::string tmp = "parsing_tmp";
    char *fileStore = new char[tmp.size()];
    strcpy(fileStore, tmp.c_str());
    store_lz77<char_type, text_offset_type>(parsing, text_length, fileStore);
    delete[] fileStore, text;
}

void hybrid::build(char* inputFilePath, uint M, char* outputFolderPath) {
    HybridSelfIndex::TRACE = false;
	HybridSelfIndex::CREATE_FMI_TEST = false;
	HybridSelfIndex::SHOW_SIZE = true;

    LZ_parsing(inputFilePath);
    std::string tmp = "parsing_tmp";
    char *fileStore = new char[tmp.size()];
    strcpy(fileStore, tmp.c_str());
    index = new HybridSelfIndex(fileStore, M, outputFolderPath);
    delete[] fileStore;

    if (std::remove(tmp.c_str()) != 0) {
        std::cerr << "Not able to delete the hybrid parsing tmp file" << std::endl;
        return;
    }
}

void hybrid::load(const std::string &inputFolderPath) {
    char *folderPath = new char[inputFolderPath.size()];
    strcpy(folderPath, inputFolderPath.c_str());
    index = new HybridSelfIndex(folderPath);
    delete[] folderPath;
}

void hybrid::locate(uchar *pat, uint m, ulong *nOcc, ulong **occ) {
    index->locate(pat, m, nOcc, occ);
}