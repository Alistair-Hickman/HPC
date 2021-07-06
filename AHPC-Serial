#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <math.h>

//Dipole spacing in Angstroms
const double dip_space = 6;
const double dip_cube = dip_space*dip_space*dip_space;
//Define Physcial constants
const double mu_0 = 4*M_PI*pow(10,-7);
const double bohr_m = 9.27*pow(10,-24);

//function for dot product of two vectors
double dot_product(std::vector<double> & vector, std::vector<double> m_hat){
 double dot = 0.0;
  for (int l=0; l<vector.size(); l++){
   dot += vector[l] * m_hat[l];
  }
 return dot;
}

int main(){

//create m vector - change depending on desired direction
std::vector<double> m_hat;
m_hat.push_back(0); //x
m_hat.push_back(0); //y
m_hat.push_back(1); //z

//Define Grid Size
const int r_x = 100;
const int r_y = 100;
const int r_z = 100;

//3 1D arrays for storing positions within grid
std::vector<double> coordinate_x;
std::vector<double> coordinate_y;
std::vector<double> coordinate_z;

//create vectors for 2D z=0 plot
std::vector<double> z0_xgrid;
std::vector<double> z0_ygrid;
std::vector<double> z0_x;
std::vector<double> z0_y;
std::vector<double> z0_value;

//Calculate max number of points required in 3D and allocate memory space
int N_x = 2*ceil(r_x/dip_space);
int N_y = 2*ceil(r_y/dip_space);
int N_z = 2*ceil(r_z/dip_space);

//reserve memory for arrays
coordinate_x.reserve(N_x*N_y*N_z);
coordinate_y.reserve(N_x*N_y*N_z);
coordinate_z.reserve(N_x*N_y*N_z);

//loop over to check if points are in ellipsoid
for (int i=0; i<N_x; i++){
 const double x = (double (i)*dip_space)/r_x-1;
  for (int j=0; j<N_y; j++){
   const double y = (double (j)*dip_space)/r_x-1;
    for (int k=0; k<N_z; k++){
      const double z = (double (k)*dip_space)/r_z-1;
      if (x*x+y*y+z*z<=1){
         coordinate_x.push_back(double(i)*dip_space);
         coordinate_y.push_back(double(j)*dip_space);
         coordinate_z.push_back(double(k)*dip_space);
      }
    }
  }
}

//store ellipsoid points in a txt file
std::ofstream ellipsoid;
ellipsoid.open("ellipsoid_points.txt");
for (unsigned int i=0; i<coordinate_x.size(); i++){
 ellipsoid << coordinate_x[i] << "\t" << coordinate_y[i] << "\t" << coordinate_z[i] << std::endl;
}
ellipsoid.close();
//print number of points to terminal as verification
std::cout << "Number of dipoles within ellipsoid=" << coordinate_x.size() << std::endl;
//create doubles to store total magnetic field components
double bx_total = 0.0;
double by_total = 0.0;
double bz_total = 0.0;
//set scalar constant
const double con = (bohr_m * mu_0)/(4*M_PI);
for (int i=0; i<coordinate_x.size(); i++){
  //set temporary holds- these build up the value for each dipole, 1 dipole per iteration
  //need to be reset to 0 at the beginning of each iteration.
   double bx_out = 0.0;
   double by_out = 0.0;
   double bz_out = 0.0;
   for (int j=0; j<coordinate_x.size(); j++){
     // generate radius components between current dipole (i) and all others (j)
     if (coordinate_x[i]!=coordinate_x[j] && coordinate_y[i]!=coordinate_y[j] && coordinate_z[i]!=coordinate_z[j]){
      //calculate radius components  between dipoles i & j
      double x_diff = coordinate_x[j] - coordinate_x[i];
      double y_diff = coordinate_y[j] - coordinate_y[i];
      double z_diff = coordinate_z[j] - coordinate_z[i];
      //calculate modulus of r
      double r_mod = sqrt((x_diff*x_diff)+(y_diff*y_diff)+(z_diff*z_diff));
      //calculate radius unit vector components
      double rx_unit = x_diff/r_mod;
      double ry_unit = y_diff/r_mod;
      double rz_unit = z_diff/r_mod;
      //store these components in a vector
      std::vector<double> r_hat;
      r_hat.push_back(rx_unit);
      r_hat.push_back(ry_unit);
      r_hat.push_back(rz_unit);
      //calculate dot product of r unit vector and m unit vector using function
      double dot = dot_product(r_hat,m_hat);
      //Apply scalar multiple of dot product to rhat and subtract mhat from result
      //Maintains a vector equation
      for (int k=0;k<r_hat.size();k++){
        r_hat[k] = ((3.0 * r_hat[k] * dot) - m_hat[k])/(r_mod*r_mod*r_mod);
      }
      //calculate components required for B_alpha
      double b_x = con*(r_hat[0]);
      double b_y = con*(r_hat[1]);
      double b_z = con*(r_hat[2]);
      //store contributions to dipole i from j in these variables.
      bx_out += b_x;
      by_out += b_y;
      bz_out += b_z;
      }
    }
 //store magnetic field for a given point by taking total contribution from all
 //other points (in b_out) and adding the points internal contribution. Temporary value required
 //for z=0 plot
 bx_total += bx_out + (con/dip_cube)*(8*M_PI*m_hat[0]/3);
 by_total += by_out + (con/dip_cube)*(8*M_PI*m_hat[1]/3);
 double bz_temp = bz_out + (con/dip_cube)*(8*M_PI*m_hat[2]/3);
 bz_total += bz_temp;

 //calculate mid point
 double z0_pos = ceil(2*r_z/dip_space) * (dip_space/2);
 //save points for z0 plot
 if (coordinate_z[i] == z0_pos){
   z0_xgrid.push_back(coordinate_x[i]);
   z0_ygrid.push_back(coordinate_y[i]);
   z0_value.push_back(bz_temp);
 }
}

// loop over saved points and verify they fall within the circle
std::vector<double> z0;
for (int m=0; m<z0_xgrid.size();m++){
 const double x_test = ((z0_xgrid[m] - r_x)*(z0_xgrid[m] - r_x));
 const double y_test = ((z0_ygrid[m] - r_y)*(z0_ygrid[m] - r_y));
 if ((x_test + y_test) <= (r_x*r_x)){
    z0_x.push_back(z0_xgrid[m]);
    z0_y.push_back(z0_ygrid[m]);
    z0.push_back(z0_value[m]);
 }
}

//create a file for the coordinates and values for z0 heat map
std::ofstream heatmap;
heatmap.open("heatmap.dat");
for (unsigned int ii=0; ii<z0_x.size(); ii++){
 heatmap << z0_x[ii] << "\t" << z0_y[ii] << "\t" << z0_value[ii] << std::endl;
}
heatmap.close();

//print number of points to terminal as verification
std::cout << "Number of dipoles in xgrid =" << z0_xgrid.size() << std::endl;
std::cout << "number in circle=" << z0_x.size() << std::endl;

//calculate average x magnetic field
double bx_avrg = (bx_total)/(coordinate_x.size()-1);
//calculate Dx
double Dx = 1-(bx_avrg/(mu_0*(bohr_m/dip_cube)));
std::cout<<"Demagnetizing factor in x=" << Dx <<std::endl;

//calculate average y magnetic field
double by_avrg = (by_total)/(coordinate_y.size()-1);
//calculate Dy
double Dy = 1-(by_avrg/(mu_0*(bohr_m/dip_cube)));
std::cout<<"Demagnetizing factor in y=" << Dy <<std::endl;

//calculate average z magnetic field
double bz_avrg = (bz_total)/(coordinate_z.size()-1);
//calculate Dz
double Dz = 1-(bz_avrg/(mu_0*(bohr_m/dip_cube)));
std::cout<<"Demagnetizing factor in z=" << Dz <<std::endl;


}
