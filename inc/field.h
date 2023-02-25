#ifndef FIELD_H
#define FIELD_H

#include "util.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <complex>

template<typename T>
class field{
    public:
        std::vector<T> data;
        int width; 
        int height;
    
        field(int h = 0, int w = 0){
            int rows = height = h;
            int cols = width = w;
            data.resize(rows*cols);
        }

        field(int h, int w, T val){
            int rows = height = h;
            int cols = width = w;

            data = std::vector<T>(rows*cols, val);
        }

        const T* operator[](const size_t index) const{
            return &data[index*width];
        }

        T* operator[](const size_t index){
            return &data[index*width];
        }

        field<T>& operator*=(const T v){
            for(int i = 0; i < width*height; i++){
                data[i] *= v;
            }
            return *this;
        }

        field<T>& operator+=(const T v){
            for(int i = 0; i < width*height; i++){
                data[i] += v;
            }
            return *this;
        }

         field<T>& operator+=(const field<T>& f){
            for(int i = 0; i < width*height; i++){
                data[i] += f.data[i];
            }
            return *this;
        }

        int v_stack(std::vector<T>& row){
            if (height == 0){
                data = row;
                height = 1;
                width = data.size();
            }else{
                if (row.size() != width) {return 1;}
                auto it = data.end();
                data.insert(it,row.begin(),row.end());
                height++;
            }
            return 0;
        }

        void fill(T v){
            for(int i = 0; i < width*height ; i++){
                data[i] = v;
            }
        }

        void msk_fill(T v, const field<int>& mask){
            for(int i = 0; i < width*height ; i++){
                if(!mask.data[i]){continue;}
                data[i] = v;
            }
        }
        
         void clamp(T min, T max){
            for(int i = 0; i < width*height ; i++){
                data[i] = std::clamp(data[i], min, max);
            }
        }

        T max() const{
            return std::max_element(data.begin(), data.end())[0];
        }

        T min() const{
            return std::min_element(data.begin(), data.end())[0];
        }

        // //save to csv
        // void save(const char* fp, char delimiter = ','){
            // std::ofstream outfile;
            // outfile.open(fp);

            // int xpos = 0;
            // for(T val : data){
                // xpos++;
                // outfile << std::setprecision(16) << val << delimiter;
                // if(xpos == width){
                // long pos = outfile.tellp();//replace last delimiter with newline
                // outfile.seekp(pos-1);
                // outfile.put('\n');
                // xpos = 0;
                // }
            // }

            // outfile.close();
        // }
        ///load from csv
        static field<T> load(const char* fp){
            auto fld = field<T>();

            std::ifstream infile;
            infile.open(fp, std::ifstream::binary);

            if(!infile.is_open()){exit(5);}
            
            for(;;){
                auto row = std::vector<T>();
                for(;;){
                    T val;
                    char delim;

                    infile >> val;
                    infile.get(delim);
                    row.push_back(val);
                    if (infile.eof()){break;}
                    if (delim == '\n' || delim == '\r'){break;}
                }
                fld.v_stack(row);
                if (infile.eof()){break;}
            }

            infile.close();
            return fld;
        }
        
        friend std::ostream& operator<< (std::ostream& os, field<T>& f){
            char delimiter = ',';
            int xpos = 0;
            for(T val : f.data){
                xpos++;
                os << val << delimiter;
                if(xpos == f.width){
                os.put('\n');
                xpos = 0;
                }
            }
            return os;
        }

        //determines if index pair is in the field
        bool in(const int i, const int j) const{
            return (0 <= i && i < height && 0 <= j && j < width);
        }

        //fold with addition
        double sum(void) const {
            double acc = 0;
            for(T val : data){
                acc += val;
            }
            return acc;
        }

        // //fold with addition excluding some elements
        // double msk_sum(const field<material>& mask) const {
            // double acc = 0;
            // for(int i = 0; i < data.size(); i++){
                // if (mask.data[i] == NO){continue;}
                // acc += data[i];
            // }
            // return acc;
        // }

        //uniform  norm
        static double norm(field<T>& a, field<T>& b);

        //subsamples a field to a lower resolution
        field<T> subsample(int level) const{
            int sample = (0x00000001 << level);

            int h = std::max(height/sample,1);
            int w = std::max(width/sample,1);

            auto subf = field<T>(h,w);

            for(int i = 0; i < h; i++){
            for(int j = 0; j < w; j++){
                subf[i][j] = (*this)[sample*i][sample*j];
            }
            }
            return subf;
        }

        //gets an element from the field clamping to be in range
        T get(int i, int j) const{
            i = std::clamp(i, 0, height-1);
            j = std::clamp(j, 0, width-1);
            return (*this)[i][j];
        }

        //gets value or returns z
        T getz(int i, int j, T z) const{
            if(0 <= i && i < height && 0 <= j && j < width){
                return (*this)[i][j];
            }else{
                return z;
            }
        }
        
        //gets value bilinearly, clamping to be in range
        T get_bilinear(double y, double x) const{
            if(0 <= y && y < height-1 && 0 <= x && x < width-1){
                double yf,xf;//fractional parts
                double yi,xi;//integral parts
                int i,j;
                
                yf = modf(y, &yi);
                xf = modf(x, &xi);
                i = (int)yi;
                j = (int)xi;
                
                T v = 0.0;
                
                v += (*this)[i  ][j  ] * (1.0 - xf) * (1.0 - yf);
                v += (*this)[i+1][j  ] * (1.0 - xf) * (      yf);
                v += (*this)[i  ][j+1] * (      xf) * (1.0 - yf);
                v += (*this)[i+1][j+1] * (      xf) * (      yf);
                
                return v;
            }else{
                int i = std::clamp(y, 0.0, height-1.0);
                int j = std::clamp(x, 0.0, width-1.0);
                return (*this)[i][j];
            }
        }

         //set value or does nothing
        void setz(int i, int j, T z){
            if(0 <= i && i < height && 0 <= j && j < width){
                (*this)[i][j] = z;
            }
        }
        
        //rotates the field //radians
        void rot(double angle){
            field<T> tmp = field<T>(height, width);
            
            double c = cos(angle);
            double s = sin(angle);
            
            for(int i = 0; i< height; i++){
            for(int j = 0; j< width; j++){
                //rotate (j,i) on to coordinate on old image (x,y)
                int jc  = j-width/2;
                int ic  = i-height/2;
                double xc = c * jc - s * ic;  
                double yc = s * jc + c * ic;  
                
                tmp[i][j] = get_bilinear(yc+height/2.0,xc+width/2.0);


            }
            }            
            data = tmp.data;
        }

        //interpolates a field from a field of a different size
        void interpolate(field<T>& f){
            int h = f.height;
            int w = f.width;

            double xscale = (double)w/(double)width;
            double yscale = (double)h/(double)height;

            for(int i = 0; i< height; i++){
            for(int j = 0; j< width; j++){

                //exact position
                double jpos = j*xscale;
                double ipos = i*yscale;

                (*this)[i][j] = get_bilinear(ipos,jpos);
            }
            }

        }
                
};

template<typename T>
int pol2cart(const field<T>& p, field<T>& c){
    double r_sf = p.width/(0.5*c.width);//radial scale factor
    for(int i = 0; i< c.height; i++){
    for(int j = 0; j< c.width; j++){
        //from rect (i,j) to polar (r,theta)
        int jc  = (j-c.width/2) ;
        int ic  = (i-c.height/2) ;
        
        int r = r_sf*sqrt(jc*jc + ic*ic);
        if(r >= p.width){continue;}

        int theta = int(p.height*(atan2(ic,jc)/(2*pi) + 1.0)) % p.height;
        
        int jw = (j+c.width/2) % c.width;
        int iw = (i+c.height/2) % c.height;
       
        c[iw][jw] = p.get_bilinear(theta,r);    
        
    }
    } 
    return 1;
}

template<typename T>
field<T> abs(const field<std::complex<T>> in){
    auto out = field<T>(in.height, in.width);
    
    for(int i = 0; i < in.data.size(); i++){
        out.data[i] = std::abs(in.data[i]);
    }
    
    return out;
}

template<typename T>
field<T> arg(const field<std::complex<T>> in){
    auto out = field<T>(in.height, in.width);
    
    for(int i = 0; i < in.data.size(); i++){
        out.data[i] = std::arg(in.data[i]);
    }
    
    return out;
}

template<typename T>
inline field<T>& operator*(field<T>& f, const T v){
    f *= v;
    return f;
}

template<typename T>
inline field<T>& operator*( const T v, field<T>& f){
    return f*v;
}

template<typename T>
inline field<T>& operator+(field<T>& f, const T v){
    f += v;
    return f;
}

template<typename T>
inline field<T>& operator+(const T v, field<T>& f){
    return f+v;
}

template<typename T>
inline  field<T>& operator+(field<T>& f, const field<T>& h){
    f += h;
    return f;
}

field<double> png2doub(const char* fp, png_mode mode = R_ONLY);
field<int> png2int(const char* fp);
void doub2png(const char* fp, const field<double>& f, png_out_mode mode = RAW);

int test_field_interpolate();

#endif FIELD_H