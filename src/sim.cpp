
#include "sim.h"

#include <random>

int project_sim(const field<double>& phantom, const scan_prop& scan, const double angle, std::default_random_engine& rng, std::vector<double>& projection){

    auto oriented = phantom;
    
    oriented.rot(angle);
    
    projection.resize(oriented.width);
    std::fill(projection.begin(), projection.end(), 0.0);
    
    for (int j = 0; j< oriented.width; j++){
        for (int i = 0; i< oriented.height; i++){
            projection[j] += oriented[i][j];
        }
        double photons;
        switch(scan.proj_int){
            case INTENSITY:
                photons = scan.I0*exp(-projection[j]*scan.att_sf);
                if(scan.noise_quanisation){
                    photons = std::round(photons);
                }
                if(scan.noise_poisson){
                    std::poisson_distribution<int> poisson(photons);
                    photons = poisson(rng);
                }
                projection[j] = photons;
                break;
            case ATTENUATION:
                break;
        }
    }
    
    
    return 1;
}



int sinogram_sim(const field<double>& phantom, const scan_prop& scan, field<double>& sinogram){
    double angle_step = 2*pi/scan.projections;

    sinogram = field<double>(scan.projections, phantom.width, 0.0);

    std::default_random_engine rng;//random engine

    auto projection = std::vector<double>(phantom.width);
    for(int i = 0; i < scan.projections; i++){
        double angle = i*angle_step;
        project_sim(phantom, scan, angle, rng, projection);
        
        std::copy(projection.begin(), projection.end(), 
                  sinogram.data.begin()+(sinogram.width*i));
        
    }
    return 1;
}

