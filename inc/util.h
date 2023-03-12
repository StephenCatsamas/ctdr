#ifndef UTIL_H
#define UTIL_H

#include<iostream>
#include<tuple>
#include<vector>
#include<random>

#include "vec.h"


extern double pi;

enum proj_integrand{
    INTENSITY,
    ATTENUATION,
};

struct scan_prop{
    proj_integrand proj_int;
    int seed;
    int I0;//photon flux (photons/pixle)
    int projections;
    double att_sf;//units for attenuation (/pixle) normally you should set this to 1/resolution to get a minium intensity of I0/e.
    bool noise_quanisation;
    bool noise_poisson;
};

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