#include <iostream>
#include <numbers>
#include <algorithm>
#include <string>
#include <complex>
#include <vector>

#include "field.h"

#include "util.h"

#include "pocketfft.h"

using namespace pocketfft;

double pi = std::numbers::pi_v<double>;

#define PROJECTIONS (512)

int project(const field<double>& phantom, const double angle, std::vector<double>& projection){
    
    auto oriented = phantom;
    
    oriented.rot(angle);
    
    projection.resize(oriented.width);
    std::fill(projection.begin(), projection.end(), 0.0);
    
    for (int j = 0; j< oriented.width; j++){
    for (int i = 0; i< oriented.height; i++){
        projection[j] += oriented[i][j];
    }
    }
    
    return 1;
}

int proj2png(char* fp, std::vector<double>& projection){
    const int projection_plot_height = 256;
    
    const int projection_size = projection.size();
    
    auto plot = field<double>(projection_plot_height, projection_size, 1);
    
    double mx = *std::max_element(projection.begin(), projection.end());
    
    for(int j = 0; j < projection_size; j++){
        int i = floor((projection_plot_height-1)*(1.0-projection[j]/mx));
        plot[i][j] = 0;
    }
    doub2png(fp, plot);
    return 1;
}

int back_project(const std::vector<double>& projection, double angle, field<double>& back_projection){
    
    for (int j = 0; j< back_projection.width; j++){
    for (int i = 0; i< back_projection.height; i++){
        back_projection[i][j] = projection[j];
    }
    } 
    back_projection.rot(-angle);
    
    return 1;
}

int recon_bp(const field<double>& phantom, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/PROJECTIONS;
    
    tomogram.fill(0.0);
    for(double angle = 0.0; angle < 2*pi; angle += angle_step){
        out.log(INF) << angle << std::endl;
        auto back_projection = field<double>(tomogram.height,tomogram.width);
        std::vector<double> projection;
        project(phantom, angle, projection);
        back_project(projection, angle, back_projection);
        tomogram += back_projection* (0.5/(PROJECTIONS*back_projection.height));
       
    }
    return 1; 
}

int recon_fbp(const field<double>& phantom, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/PROJECTIONS;   
    
    tomogram.fill(0.0);
    for(double angle = 0.0; angle < 2.0*pi; angle += angle_step){
        out.log(INF) << angle << std::endl;
        auto back_projection = field<double>(tomogram.height,tomogram.width);
        std::vector<double> projection;
        project(phantom, angle, projection);
        
        //fft stuff
        //pad projection
        projection.resize(2*projection.size(), 0.0);

        shape_t shape = {projection.size()};
        stride_t stride_r = {sizeof(double)};
        stride_t stride_c = {sizeof(std::complex<double>)};
        shape_t axes = {0};

        auto f_projection = std::vector<std::complex<double>>(shape[0]/2 + 1);
        double scale_factor = 0.5/projection.size();

        r2c(shape,
            stride_r, 
            stride_c, 
            axes,
            FORWARD, 
            projection.data(), 
            f_projection.data(), 
            scale_factor);
        
        //filter
        for(int i = 0; i < shape[0]/2 + 1; i++){
            double weight = 2*pi*(i);
            // double window = 0.54+0.46*cos(2*pi*i/(4*shape[0]/2));
            f_projection[i] *= weight;
        }    
        
        c2r(shape,
            stride_c, 
            stride_r, 
            axes,
            BACKWARD, 
            f_projection.data(), 
            projection.data(), 
            1.0);
   
        //unpad projection
        projection.resize(back_projection.width);
        //end fft stuff    
                
        back_project(projection, angle, back_projection);
        tomogram += back_projection * (0.5/(PROJECTIONS*back_projection.height));
       
    }
    return 1; 
}


int recon_dfi(const field<double>& phantom, field<double>& tomogram){
    
    const double angle_step = 2.0*pi/PROJECTIONS;   
    auto f_polar_proj = field<std::complex<double>>(PROJECTIONS,tomogram.width + 1);
    auto f_tomogram = field<std::complex<double>>(tomogram.height,tomogram.width);
    tomogram.fill(0.0);
    for(int i = 0; i < PROJECTIONS; i++){
        double angle = i * angle_step;
        // out.log(INF) << angle << std::endl;
        
        std::vector<double> projection;
        project(phantom, angle, projection);
        
        projection.resize(2*projection.size(), 0.0);
        shape_t shape = {projection.size()};
        stride_t stride_r = {sizeof(double)};
        stride_t stride_c = {sizeof(std::complex<double>)};
        shape_t axes = {0};

        auto f_projection = std::vector<std::complex<double>>(shape[0]/2 + 1);
        double scale_factor = 0.5/projection.size();

        r2c(shape,
            stride_r, 
            stride_c, 
            axes,
            FORWARD, 
            projection.data(), 
            f_projection.data(), 
            scale_factor);
        
        
        std::copy(f_projection.begin(), f_projection.end(), 
                  f_polar_proj.data.begin()+(f_polar_proj.width*i));       
    }
    out.log(INF) << "interpolaring to cart" << std::endl;
    
    //interpolar polar to rect field
    int height = f_tomogram.height;
    int width = f_tomogram.width;
    for(int i = 0; i< height; i++){
    for(int j = 0; j< width; j++){
        //from rect (i,j) to polar (r,theta)
        int jc  = (j-width/2) ;
        int ic  = (i-height/2) ;
        
        double r = 2.0*sqrt(jc*jc + ic*ic);
        double dump;
        double theta = PROJECTIONS*modf(atan2(ic,jc)/(2*pi) + 1.5, &dump);
        
        //exact position
        double js;
        double is;
        double jfrac = modf(r, &js);
        double ifrac = modf(theta*PROJECTIONS, &is);
        // std::cout << << is << ifrac << js << jfrac << std::endl;

        //bilinear interpolation
        std::complex<double> a = (ifrac)    *(1.0-jfrac)*f_polar_proj.get(is+1,js) ;
        std::complex<double> b = (ifrac)    *(jfrac)    *f_polar_proj.get(is+1,js+1) ;
        std::complex<double> c = (1.0-ifrac)*(1.0-jfrac)*f_polar_proj.get(is,js) ;
        std::complex<double> d = (1.0-ifrac)*(jfrac)*f_polar_proj.get(is,js+1);
        

        int jw  = (j+width/2 ) % width;
        int iw  = (i+height/2 ) % height;
        f_tomogram[iw][jw] = a+b+c+d;   

        // f_tomogram[i][j] = phantom[i][j];     
        
    }
    } 
    out.log(INF) << "fourier transforming" << std::endl;
    
    shape_t shape_2d = {(size_t)f_tomogram.width, (size_t)f_tomogram.height};
    stride_t stride_2d = {sizeof(std::complex<double>),sizeof(std::complex<double>)*f_tomogram.width};
    shape_t axes_2d = {0,1};   
               
    c2c(shape_2d,
        stride_2d, 
        stride_2d, 
        axes_2d,
        BACKWARD, 
        f_tomogram.data.data(), 
        f_tomogram.data.data(), 
        1.0);
        
    auto f_polar_proj_r = arg(f_polar_proj);
    out.log(INF) << "saving sinogram" << std::endl;
    double mx = f_polar_proj_r.max();
    double mn = f_polar_proj_r.min();
    f_polar_proj_r += -mn;
    f_polar_proj_r *= 1/(mx-mn);
    out.log(INF) << mx << std::endl;
    doub2png("data/sinogram.png", f_polar_proj_r);
    auto f_tomogram_r = abs(f_tomogram);
    out.log(INF) << "saving f_tomogram" << std::endl;
    mx = f_tomogram_r.max();
    mn = f_tomogram_r.min();
    f_tomogram_r += -mn;
    f_tomogram_r *= 1/(mx-mn);
    out.log(INF) << mx << std::endl;
    doub2png("data/f_tomogram.png", f_tomogram_r);
    out.log(INF) << "saving tomogram" << std::endl;
    mx = tomogram.max();
    tomogram *= 1/mx;
    out.log(INF) << mx << std::endl;
    doub2png("data/tomogram.png", tomogram);
 
    return 1;
}


int fft_2d_test(){
    auto r = field<double>(128,128, 0.0);
    auto f = field<std::complex<double>>(
                        r.height/2 + 1,
                        r.width, 0.0);
    
    f[10][10] = 1.0;
    
    shape_t shape = {(size_t)r.width, 
                     (size_t)r.height};
    stride_t stride_r = {sizeof(double),
                         sizeof(double)*r.width};
    stride_t stride_c = {sizeof(std::complex<double>), 
                         sizeof(std::complex<double>)*r.width};
    shape_t axes = {1};   
               
    c2r(shape,
        stride_c, 
        stride_r, 
        axes,
        BACKWARD, 
        f.data.data(), 
        r.data.data(), 
        1.0);

    auto f_abs = abs(f);
    out.log(INF) << "saving fourier" << std::endl;
    double mx = f_abs.max();
    f_abs *= 1/mx;
    out.log(INF) << mx << std::endl;
    doub2png("data/fourier.png", f_abs);
    
    out.log(INF) << "saving real" << std::endl;
    mx = r.max();
    double mn = r.min();
    r += -mn;
    r *= 1/(mx - mn);
    out.log(INF) << mx << std::endl;
    doub2png("data/real.png", r);
    return 1;
}

int main() {
    out.log(INF) << "ctdr" << std::endl;
    
    // fft_2d_test();
    
    out.log(INF) << "loading phantom" << std::endl;
    
    auto phantom = png2doub("data/phantom.png", R_ONLY);
    
    auto tomogram = field<double>(phantom.height, phantom.width);
        
    // recon_bp(phantom, tomogram);
    // recon_fbp(phantom, tomogram);
    recon_dfi(phantom, tomogram);
    
    double tm_mx = tomogram.max();
    out.log(INF) << tm_mx << std::endl;
    tomogram.clamp(0.0,1.0);
    
    out.log(INF) << "saving tomogram" << std::endl;
    doub2png("data/tomogram.png", tomogram);

    return 0;
}
