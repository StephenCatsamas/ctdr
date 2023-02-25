#ifndef UTIL_H
#define UTIL_H

#include<iostream>
#include<tuple>
#include<vector>

enum png_mode{
    R_ONLY,
    RGBA,
};

enum png_out_mode{
    RAW,
    SCALE,
    SCALE_ZERO,
};

enum print_mode{
    DBG,
    INF,
    WRN,
    ERR,
};

class logger{
    public:
        print_mode level;

        std::ostream* printer;

        logger(print_mode, std::ostream& = std::cout);

        std::ostream& log(print_mode);
};

extern logger out;

std::tuple<std::vector<unsigned char>,unsigned,unsigned> decode_png(const char* filename);
void encode_png(const char* , std::vector<unsigned char>& , unsigned , unsigned);

#endif UTIL_H