#ifndef SIM_H
#define SIM_H

#include "util.h"
#include "field.h"


int sinogram_intensity(const field<double>& phantom, int projections, field<double>& sinogram);
int sinogram_attenuation(const field<double>& phantom, int projections, field<double>& sinogram);

int project_intensity(const field<double>& phantom, const double angle, std::vector<double>& projection);
int project_attenuation(const field<double>& phantom, const double angle, std::vector<double>& projection);



#endif SIM_H