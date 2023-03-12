#ifndef RECON_H
#define RECON_H

#include "util.h"
#include "field.h"

int project(const field<double>& phantom, const double angle, std::vector<double>& projection);
int back_project(const std::vector<double>& projection, double angle, field<double>& back_projection);


int recon_bp(const field<double>& sinogram, const int n_proj, field<double>& tomogram);
int recon_fbp(const field<double>& sinogram, const int n_proj, field<double>& tomogram);
int recon_dfi(const field<double>& sinogram, const int n_proj, const int pf, field<double>& tomogram);
int recon_art(const field<double>& sinogram, const int n_proj, field<double>& tomogram);

int proj2png(char* fp, std::vector<double>& projection);



#endif RECON_H