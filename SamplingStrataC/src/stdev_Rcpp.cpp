/* standard deviation calculated with distances */

# include <Rcpp.h>
# include <iostream>
//# include <vector>
using std::cout;
using std::cin;
using std::endl;
//using std::vector;

using namespace Rcpp;


// [[Rcpp::export]]
double stdev_Rcpp(DataFrame dataset, 
              int i,
              double fitting,
              double range,
              int kappa) {
  NumericVector LON = dataset["LON"];
//  cout << "LON: " << LON << endl;
  NumericVector LAT = dataset["LAT"];
//  cout << "LAT: " << LAT << endl;

  std::ostringstream str1; 
  str1 << i; 
  std::string str2 = str1.str();
  
  std::string s1("Y");
  std::string y = s1 + str2;
  std::string s2("VAR");
  std::string v = s2 + str2;
  // std::cout << "y: " << y << std::endl;
  // std::cout << "v: " << v << std::endl;
  
  /*
  String y("Y");
  y += "2";
  String v("VAR");
  v += "2";
 */

  NumericVector target = dataset[y];
  // std::cout << "Y1: " << target << std::endl;
  NumericVector var = dataset[v];
  // cout << "VAR1: " << var << endl;
  int elems = LON.length();
  // cout << "elems: " << elems << endl;
  NumericVector dist(elems*elems);
  NumericVector z_z(elems*elems);
  NumericVector somma_coppie_var(elems*elems);
  NumericVector prod_coppie_std(elems*elems);
  NumericVector spatial_cov(elems*elems);
  double sd1 = 0;
  double sd2 = 0;
  double sum_z_z = 0;
  double sum_somma_coppie_var = 0;
  double var_strato = 0;
  double stdev_strato;

  for ( int i=0; i<elems; i++ ) {
    for ( int j=0; j<elems; j++ ) {
      dist[i*elems+j] = sqrt((LON[i]-LON[j])*(LON[i]-LON[j]) + (LAT[i]-LAT[j])*(LAT[i]-LAT[j]));
      z_z[i*elems+j] = (target[i] - target[j]) * (target[i] - target[j]);
      somma_coppie_var[i*elems+j] = (var[i] + var[j]);
      prod_coppie_std[i*elems+j] = sqrt(var[i] * var[j]);
      //    cout << "distance: " << dist[i*elems+j] << endl;
      //    cout << "z_z: " << z_z[i*elems+j] << endl;
      //    cout << "somma_coppie_var: " << somma_coppie_var[i*elems+j] << endl;
      //    cout << "prod_coppie_std: " << prod_coppie_std[i*elems+j] << endl;
    }
  }
  for ( int i=0; i<elems; i++ ) {
    for ( int j=0; j<elems; j++ ) {
      spatial_cov[i*elems+j] = prod_coppie_std[i*elems+j]*exp(-kappa*dist[i*elems+j]/range);
      //    cout << "spatial_cov: " << spatial_cov[i*elems+j] << endl;
    }
  }
for ( int i=0; i<elems; i++ ) {
  for ( int j=0; j<elems; j++ ) {
      sum_z_z = sum_z_z + z_z[i*elems+j];
      sum_somma_coppie_var = sum_somma_coppie_var + somma_coppie_var[i*elems+j] -2 * spatial_cov[i*elems+j];
    }
  }
  //    cout << "sum_z_z: " << sum_z_z << endl;
  //    cout << "sum_somma_coppie_var: " << sum_somma_coppie_var << endl;
  if (elems == 1) {
    //somma_coppie_var[1] = 0;
    //spatial_cov[1] = 0;
    var_strato = 0;
  }
  if (elems > 1) {
    sd1 = sqrt(sum_z_z / (2*elems*elems));
    sd2 = sqrt(sum_somma_coppie_var / (2*elems*elems));
    var_strato = (sd1*sd1 / fitting) + sd2*sd2;
    //    cout << "sd1: " << sd1 << endl;
    //    cout << "sd2: " << sd2 << endl;
    //    cout << "var_strato: " << var_strato << endl;
  } 
  stdev_strato = sqrt(var_strato);
  return stdev_strato;
}

