#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(Rcpp)]]
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
#include <algorithm>
#include <string>
#include <iostream>
#include <math.h>
using namespace Rcpp;
using namespace std;
using namespace arma;
using namespace RcppParallel;


//colnames(stratif) <- toupper(colnames(stratif))
//  colnames(errors) <- toupper(colnames(errors))
//  checkData(strata = stratif, errors = errors)
//  ordina_variabili <- function(dati, prefisso, n_var) {
//    if (!is.data.frame(dati)) 
//      stop()
//     as.matrix(dati[, paste(prefisso, 1:n_var, sep = ""), 
//                drop = FALSE])
//  }

//ordina_variabili <- function(dati, prefisso, n_var) {
//  if (!is.data.frame(dati)) 
//    stop()
//    as.matrix(dati[, paste(prefisso, 1:n_var, sep = ""), 
//              drop = FALSE])
//}

// [[Rcpp::export]]
CharacterVector uppercasevec(DataFrame df) {
  vector<string> strings = df.attr("names");
  int len = strings.size();
  
  for( int i=0; i < len; i++ ) {
    std::transform(strings[i].begin(), strings[i].end(), strings[i].begin(), ::toupper);
  }
  
  return wrap(strings);
  
}

// [[Rcpp::export]]
int  grepRcpp(const std::string &pattern, StringVector &x ) {
  std::string parentstring;
  int len = x.size();
  int total=0;
  for(int i = 0; i < len; ++i) {
    parentstring = x[i];
    if (parentstring.find(pattern)){}else{total+=1;}
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::IntegerVector seqRcpp(int vx){
  Rcpp::IntegerVector invec(vx);
  Rcpp::IntegerVector outvec(vx);
  for(int i=0; i<invec.size(); i++) {
    outvec[i] = i+1;
  }
  return outvec;
}
// [[Rcpp::export]]
std::string int_to_string(int &addr){
  stringstream ss;
  ss << addr;
  std::string str = ss.str();
  return str;
}
//thanks to https://stackoverflow.com/questions/43182003/concatenate-stringvector-with-rcpp
// [[Rcpp::export]]
CharacterVector paste3(CharacterVector lhs,CharacterVector rhs)
{
  using proxy_t = Rcpp::internal::string_proxy<STRSXP>;
  
  std::vector<std::string> res(lhs.begin(), lhs.end());
  std::transform(res.begin(), res.end(), rhs.begin(), res.begin(),
                 [&](const std::string& x, const proxy_t& y) {
                   return x + y;
                 }
  );
  
  return  wrap(res) ;
}
//thanks to https://stackoverflow.com/questions/28442582/reproducing-r-rep-with-the-times-argument-in-c-and-rcpp
struct Sum : Worker {
  const RVector<int> input;
  int value;
  Sum(const IntegerVector& input) : input(input), value(0) {}
  Sum(const Sum& sum, Split) : input(sum.input), value(0) {}
  void operator()(std::size_t begin, std::size_t end) {
    value += std::accumulate(input.begin() + begin, input.begin() + end, 0);
  }
  void join(const Sum& rhs) {
    value += rhs.value;
  }
};

struct Fill: Worker {
  const RVector<double> input;
  const RVector<int> times;
  RVector<double> output;
  std::size_t ind;
  Fill(const NumericVector& input, const IntegerVector& times, NumericVector& output)
    : input(input), times(times), output(output), ind(0) {}
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ind += times[i], ++i)
      std::fill(output.begin() + ind, output.begin() + ind + times[i], input[i]);
  }
};

// [[Rcpp::export]]
NumericVector repRcpp(NumericVector x,  IntegerVector times) {
  std::size_t n = times.size();
  if (n != 1 && n != x.size())
    stop("Invalid 'times' value");
  Sum s(times);
  parallelReduce(0, n, s);
  NumericVector res = no_init(s.value);
  Fill f(x, times, res);
  parallelFor(0, n, f);
  return res;
}

// [[Rcpp::export]]
NumericVector repRcpp2(int x, int n) {
  NumericVector myvector(n);
  int ind = 0;
  for (int i = 0; i < n; ++i) {
    myvector[ind] = x;
    ind = ind + 1;
  }
  return myvector;
}

// [[Rcpp::export]]
NumericVector repRcppdouble(double x, int n) {
  NumericVector myvector(n);
  int ind = 0;
  for (int i = 0; i < n; ++i) {
    myvector[ind] = x;
    ind = ind + 1;
  }
  return myvector;
}
//[[Rcpp::export]]
NumericMatrix testDFtoNM1(DataFrame x) {
  int nRows=x.nrows();  
  NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y(_,i)=NumericVector(x[i]);
  }  
  return y;
}

// [[Rcpp::export]]
NumericMatrix ordina_variabiliRcpp(DataFrame &dati,CharacterVector prefisso,int &n_var) {
  std::string y =  dati.attr("class");
  std::string z = "data.frame";
  if (y!=z){stop("dati is not a data set");};
  CharacterVector nvar;
  int number;
  CharacterVector subset(n_var);
  StringVector  output(n_var);
  for(int i = 0; i < n_var; ++i) {
    number = i+1;
    nvar = int_to_string(number);
    subset= paste3(prefisso,nvar);
    output[i]=as<string>(subset);
  }  
  DataFrame outdf = dati[output];
  NumericMatrix outmat = testDFtoNM1(outdf);
  return outmat ;
}
// [[Rcpp::export]]
StringVector ndomRcpp(CharacterVector prefisso,int n_vars) {
  CharacterVector nvar;
  int number;
  CharacterVector subset(n_vars);
  StringVector  output(n_vars);
  for(int i = 0; i < n_vars; ++i) {
    number = i+1;
    nvar = int_to_string(number);
    subset= paste3(prefisso,nvar);
    output[i]=as<string>(subset);
  }  
  
  return output ;
}
// [[Rcpp::export]]
IntegerVector nvaluesRcpp(StringVector nom_dom, NumericMatrix dom ){
  IntegerVector levs;
  for(int i = 0; i < nom_dom.size(); ++i) {
    levs = sort_unique(dom(_,i));
  }  
  //levs.attr("names") = nom_dom;
  return levs.size();
}

// [[Rcpp::export]]
arma::mat crea_disjRcpp(DataFrame data, CharacterVector pref){
  
  DataFrame col= data[pref];
  NumericMatrix col1 = testDFtoNM1(col);
  NumericVector levs= sort_unique(col1(_,0));
  arma::mat out(data.nrow(),levs.size());
  for (int j=0; j<levs.size();j++) {
    for (int i=0; i<data.nrow();i++) {
      if(col1(i,0)==levs[j])
        out(i,j)=1;
      else{out(i,j)=0;};
    }
  }
  return out;
  
}
// [[Rcpp::export]]
CharacterVector dfrow(DataFrame& x, int num) {
  // Suppose I need to get the 10th row
  int nCols=x.size();
  CharacterVector y(nCols);
  for (int j=0; j<nCols;j++) {
    CharacterVector column = x[j] ;
    y[j] = column[num] ;
  }
  return y;
}  
// [[Rcpp::export]]
CharacterVector int_to_charvec(int addr){
  stringstream ss;
  ss << addr;
  std::string str = ss.str();
  return wrap(str);
}
// [[Rcpp::export]]
std::string paste4(CharacterVector lhs,CharacterVector rhs)
{
  using proxy_t = Rcpp::internal::string_proxy<STRSXP>;
  
  std::vector<std::string> res(lhs.begin(), lhs.end());
  std::transform(res.begin(), res.end(), rhs.begin(), res.begin(),
                 [&](const std::string& x, const proxy_t& y) {
                   return x + y;
                 }
  );
  
  return Rcpp::as<string>(wrap(res));
}
// [[Rcpp::export]]
NumericMatrix testDFtoNM3(DataFrame x) {
  int nRows=x.nrows();  
  NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y(_,i)=NumericVector(x[i]);
  }  
  //y.attr("names")=x.attr("names");
  return wrap(y);
}

// [[Rcpp::export]]
NumericMatrix ordina_variabiliRcpp2(DataFrame dati,CharacterVector prefisso,int n_var) {
  CharacterVector nvar(n_var);
  int number=0;
  std::string subset;
  StringVector    output(n_var);
  for(int i = 0; i < n_var; ++i) {
    number = i+1;
    nvar = int_to_charvec(number);
    subset= paste4(prefisso,nvar);
    output[i]=subset;
  }  
  DataFrame outdf = dati[output];
  NumericMatrix outmat = testDFtoNM3(outdf);
  return  outmat;
}

// [[Rcpp::export]]
arma::mat cbindRcpp(arma::mat a,arma::mat b) {
  arma::mat out = join_rows(a, b);
  return out ;
}

// [[Rcpp::export]]
NumericMatrix cbindRcpp1(NumericMatrix a, NumericMatrix b) {
  int acoln = a.ncol();
  int bcoln = b.ncol();
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b(_, j - acoln);
    }
  }
  return out;
}
// [[Rcpp::export]]
arma::rowvec vecconcat(arma::rowvec x,arma::rowvec y) {
  arma::rowvec z(x.size() + y.size());
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  return(z);
}                      

// [[Rcpp::export]]
CharacterVector stringrepRcpp(StringVector x, int y) {
  std::vector<std::string> myvector(y);
  int ind = 0;
  for (int j = 0; j < y; ++j) {
    myvector[ind] = x[0];
    ind = ind + 1;
  }
  return wrap(myvector);
}
// [[Rcpp::export]]
std::vector<std::string> concat(std::vector<std::string> x,
                                std::vector<std::string> y) {
  
  std::vector<std::string> z(x.size() + y.size());
  
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  
  return(z);
}

// [[Rcpp::export]]
std::vector<std::string>  stringrepRcpp2(StringVector x, int y) {
  std::vector<std::string> myvector(y);
  int ind = 0;
  for (int j = 0; j < y; ++j) {
    myvector[ind] = x[0];
    ind = ind + 1;
  }
  return myvector;
}

// [[Rcpp::export]]
arma::mat crea_aRcpp(arma::vec  N11, arma::mat s11, 
                     arma::vec nocens11, arma::mat  m11,
                     arma::vec cv11, double epsilon= 1e-11){ 
  int n;
  arma::mat  v=s11;
  arma::mat s2(s11);
  int sn = s11.n_rows;int sc=s11.n_cols;
  for (int i = 0; i < sn; i++) {
    for (int j = 0; j < sc; j++) {
      s2(i,j)=s11(i,j) * v(i,j) ;
    }
  }
  n = N11.size();
  arma::vec N2(n);
  for(int j = 0; j < n; j++){
    N2[j]=N11[j]*N11[j];
  }
  //numA <- (N^2) * (s^2) * nocens
  arma::vec N2nocens(n);
  for(int j = 0; j < n; j++){
    N2nocens[j]=N2[j]* nocens11[j];
  }
  arma::mat numA(s2);
  int s2c = s2.n_cols;
  
  for(int i = 0; i < s2c; i++){
    for(int j = 0; j < n; j++){
      numA(j,i)= s2(j,i)*N2nocens[j];
    }}//ok 27/12/2017
  
  mat Nmcv= (m11.t()*N11)*(cv11.t());//ok 
  mat denA1 = pow(Nmcv.col(0),2); 
  //denA2 <- colSums(N * (s^2) * nocens)
  arma::vec N3(n);
  for(int j = 0; j < n; j++){
    N3[j]=N11[j]* nocens11[j];
  }//correct
  
  arma::mat denA2=s2.t()*N3;//ok
  //  denA <- denA1 + denA2 + epsilon
  arma::vec denA=denA1 + denA2+ epsilon; //ok
  
  int mc =numA.n_cols;
  n = numA.n_rows;
  arma::mat a11(n,mc);
  
  for(int j = 0; j < n; j++){
    for(int i = 0; i < mc ; i++)
    {
      a11(j,i)= numA(j,i)/denA[i];
    }
  } 
  return a11;
}
// [[Rcpp::export]]
arma::mat a13 (arma::mat x) {
  return(x.cols(1, 1)) ;
}
// [[Rcpp::export]]
NumericVector chromyNew1(double& alfatot, double& diff, int& iter,
                         arma::vec alfa, arma::vec alfanext, 
                         arma::vec x,const arma::mat a,const 
                         arma::vec cost, const double& epsilon = 1e-11,
                         const int& maxiter=200,const  bool& realAllocation=true){
  arma::vec den1; 
  double den2; 
  int ncols =a.n_cols;
  int nrows = a.n_rows;
  int nRows;
  arma::mat out(nrows,ncols);
  arma::mat tout(ncols,nrows);
  arma::vec sqrt_rsums;   arma::mat at(a.t()); arma::mat atout(ncols,nrows); arma::mat t_atout(a.t());
  arma::vec alfa2; 
  nRows = a.n_rows;
  arma::vec out2(nRows);
  arma::vec sqrtcost=sqrt(cost);
  arma::vec den1den2eps = (den1*den2)+epsilon;
  arma::mat  ax = a.t() * x;
  arma::mat   v=ax;
  arma::mat ax2(ax);
  
  int sn = ax.n_rows;
  int sc=ax.n_cols;
  arma::mat alfaax2(sn,sc);
  arma::vec out3(ncols); 
  //  mat alfaax2(nrows1,ncols1);
  arma::mat ax2alfatot(ax2);
  arma::mat alfanextmat(sn,sc);
  NumericMatrix alfanextM; 
  NumericVector alfanextV; 
  while (diff > epsilon && iter < maxiter) {
    
    iter = iter + 1;
    den1 = sqrt(a*alfa);//ok 28/12/2017
    for(int j = 0; j < nrows; j++){
      for(int i = 0; i < ncols ; i++)
      {
        out(j,i)= a(j,i)*cost[j];
      }
    } 
    tout=out.t();
    //remember matrix has been transposed so num rows and num cols go in opposite places
    for(int j = 0; j <ncols; j++){
      for(int i = 0; i < nrows ; i++)
      {
        atout(j,i)= tout(j,i)*alfa[j];
      }
    } 
    t_atout= atout.t();
    
    
    for(int i = 0; i < nRows; i++){
      out2(i) = sum(t_atout.row(i));
    }
    
    sqrt_rsums =sqrt(out2); 
    den2 = sum(sqrt_rsums);
    sqrtcost=sqrt(cost);
    den1den2eps = (den1*den2)+epsilon;
    x = sqrtcost/den1den2eps;
    ax = a.t() * x;
    v=ax;
    
    for (int i = 0; i < sn ; i++) {
      for (int j = 0; j < sc ; j++) {
        ax2(i,j)=ax(i,j) * v(i,j) ;
      }
    }
    
    for(int j = 0; j < sn; j++){
      for(int i = 0; i < sc ; i++)
      {
        alfaax2(j,i)= ax2(j,i)*alfa[j];
      }
    } 
    
    for(int i = 0; i < ncols; i++){
      out3(i) = sum(alfaax2.row(i));
    }
    alfatot=sum(out3);
    if(alfatot == 0){alfatot= epsilon;}
    
    ax2alfatot =ax2/alfatot;
    
    
    for(int j = 0; j < sn; j++){
      for(int i = 0; i < sc ; i++)
      {
        alfanextmat(j,i)= ax2alfatot(j,i)*alfa[j];
      }
    } 
    alfanextM = wrap(alfanextmat);
    alfanextV = alfanextM(_,0);
    alfanext= as<vec>(alfanextV);
    diff = max(abs(alfanext - alfa));
    alfa = alfanext;
    alfa2 = alfanext;
    
  }
  
  arma::vec n;
  
  if (realAllocation == false) {n =ceil(1/x);}
  else if (realAllocation == true) {n =1/x;};
  return  wrap(n);
}


// [[Rcpp::export]]
arma::vec bethelRcppOpen(DataFrame &stratif, 
                         DataFrame &errors,
                         int minnumstrat = 2, 
                   int  maxiter = 200, 
                   int maxiter1 = 25, 
                   LogicalVector printa=true,
                   LogicalVector realAllocation=true, 
                   double epsilon = 1e-11){
  
  CharacterVector strat = uppercasevec(stratif);
  stratif.attr("names")= strat;
  CharacterVector err = uppercasevec(errors);
  errors.attr("names") = err;
  
  int iter1 =0;
  //  val <- NULL
  arma::mat val;
  
  arma::mat s;
  arma::rowvec cv;
  int  nstrat = stratif.nrow();
  int    nvar = grepRcpp("CV",err);
  int    ndom = errors.nrow();
  Rcpp::IntegerVector varloop= seqRcpp(nvar);
  Rcpp::IntegerVector strloop= seqRcpp(nstrat);
  NumericMatrix med = ordina_variabiliRcpp(stratif, "M", nvar);
  NumericMatrix  esse = ordina_variabiliRcpp(stratif, "S", nvar);
  if (med.ncol() != esse.ncol())
    stop("Error: Number of M variables don't match the number of S variables");
  if (med.ncol() != nvar) 
    stop("Error: Number of variables don't match the number of planned CV");
  arma::rowvec N =  Rcpp::as<arma::rowvec>(stratif["N"]);
  arma::rowvec cens =  Rcpp::as<arma::rowvec>(stratif["CENS"]);
  

  for(int i = 0;i <  N.size(); ++i){
    if (N[i] < minnumstrat) cens[i] = 1;
  }
  arma::vec  cost = Rcpp::as<vec>(stratif["COST"]);
  if (cost.empty()){cost = repRcpp(1, nstrat);};
  
  if (cens.empty()){cens = repRcpp(0, nstrat);};
  arma::rowvec nocens = 1 - cens;
  StringVector nom_dom = ndomRcpp("DOM",ndom);
  NumericMatrix dom = ordina_variabiliRcpp(stratif, "DOM", ndom);
  IntegerVector nvalues =  nvaluesRcpp(nom_dom,dom );
  
  CharacterVector cvDom;
  
  
  Rcpp::IntegerVector ndomvalues;
  NumericMatrix err_row;
  
  NumericMatrix outmat;
  outmat = ordina_variabiliRcpp2(errors, "CV", nvar);
  arma::mat armaout = as<mat>(outmat);
  arma::rowvec cvx;
  int nvaluesint=nvalues.size();
  std::vector<std::string> cvDom2;
  std::string cvd2;
  std::vector<std::string> cvdom(nvar);
  int k2;
  for (int k=0; k<ndom;k++) {
    ndomvalues =seqRcpp(nvaluesint);
    cvx=outmat(k,_);
    for (int k1=0; k1<nvaluesint;k1++) {
      cv = vecconcat(cv,cvx);
    }
  }
  
  //for (i in 1:nc) {
  //  m <- cbind(m, disj[, i] * med)
  //  s <- cbind(s, disj[, i] * esse)
  //}
  arma::mat disj = crea_disjRcpp(stratif,nom_dom);
  int nc = disj.n_cols;
  arma::mat medmat = as<mat>(med);
  arma::mat essemat = as<mat>(esse);
  arma::mat mprep;
  int mc =medmat.n_cols;
  int mr = medmat.n_rows;
  arma::mat m;
  arma::vec disjvec;
  arma::mat mout(mr,mc);
  arma::mat sout(mr,mc);
  //for (i in 1:nc) {
  //  m <- cbind(m, disj[, i] * med)
  //  s <- cbind(s, disj[, i] * esse)
  //}
  for (int i=0; i< nc ;i++) {
    disjvec=disj.cols(i, i);
    for(int j = 0; j < mc; j++){
      for(int k = 0; k < mr ; k++)
      {
        mout(k,j)= (medmat(k,j)*disjvec[k]);
        sout(k,j)= (essemat(k,j)*disjvec[k]);
      }
    }
    m= join_rows(m, mout);
    s= join_rows(s, sout);
  }
  cvDom=stringrepRcpp(nom_dom,cvx.size()*nvaluesint);
  for (int k1=0; k1<nvaluesint;k1++) {
    k2=k1+1;
    cvd2=int_to_string(k2);
    cvdom=stringrepRcpp2(cvd2,nvar);
    cvDom2=concat(cvDom2,cvdom);
  }
  int nvar1=cv.size();
  double nvar2 =cv.size();
  varloop = seqRcpp(nvar1);
  NumericVector varfin=repRcpp2(0,nvar1);
  NumericVector totm = repRcpp2(0, nvar1);
  NumericVector alfa2  = repRcpp2(0, nvar1);
  arma::mat a1(s.n_rows,s.n_cols);
  arma::vec cvvec =cv.t();
  arma::vec Nvec= N.t();
  arma::vec nocensvec=nocens.t();
  a1 = crea_aRcpp(Nvec, s, nocensvec, m,  cvvec, epsilon);
  double alfatot1=0;
  double diff1=999;
  int iter= 0;
  double oneovernvar= 1/nvar2;
  arma::vec alfa1=repRcppdouble(oneovernvar, nvar1);
  arma::vec alfanext1  = repRcpp2(0, nvar1);
  arma::vec x1= repRcppdouble(0.1, nstrat);
  
  //n <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, nvar)), array(0.1, dim = c(nstrat, 1)))
  //chromyNew2(double& alfatot, double& diff, int& iter,
  //vec& alfa, vec& alfanext, vec& x, mat& a,vec& cost, const double& epsilon = 1e-11,
  //  const int& maxiter=200,const bool& realAllocation=false)
  arma::vec n = chromyNew1(alfatot1,diff1,iter, alfa1, alfanext1,x1,a1,cost,epsilon,
                     maxiter,realAllocation);
  
  int   contx = sum(n > Nvec);
  
  for(int i = 0;i <  N.size(); ++i){
    cens[i] = 0;
    if (n[i] > Nvec[i]) cens[i] = 1;
  }
  nocens = 1 - cens;
  
  
  arma::vec outvec(N.size());
  for(int i = 0;i <  N.size(); ++i){
    if(minnumstrat < N[i]){outvec[i]= minnumstrat;}else{outvec[i]=N[i];};
  }
  
  for(int i = 0;i <  N.size(); ++i){
    if (n[i] < minnumstrat) {n[i]=outvec[i];};
  }
  //cout << "outvec: " << outvec << "\n  ";
  
  while (contx > 0 && iter1 < maxiter1) {
    iter1 = iter1 + 1;
    //a <- crea_a()
    a1 = crea_aRcpp(Nvec, s, nocensvec, m,  cvvec, epsilon);
    //n <- chromy(0, 999, 0, c(rep(1/nvar, nvar)), c(rep(0, nvar)), array(0.1, dim = c(nstrat, 1)))
    arma::vec n = chromyNew1(alfatot1,diff1,iter, alfa1, alfanext1,x1,a1,cost,epsilon,
                   maxiter,realAllocation);;

 
    
    //contx <- sum(n > N)
    contx = sum(n > Nvec);
    //cens[n > N] <- 1
    //cout << "Nvec: " << Nvec << "\n";
    for(int i = 0;i <  N.size(); ++i){
      cens[i] = 0;
      if (n[i] > Nvec[i]) cens[i] = 1;
      //cout << "\n n: " << n[i] << "  Nvec: " << Nvec[i] << "  cens: " << cens[i];
    }
    
    nocens = 1 - cens;
    // n <- check_n()
    for(int i = 0;i <  N.size(); ++i){
      if(minnumstrat < Nvec[i]){outvec[i]= minnumstrat;}else{outvec[i]=Nvec[i];};
    }
    
    for(int i = 0;i <  N.size(); ++i){
      if (n[i] < minnumstrat) {n[i]=outvec[i];};
    }
  }
  
  
  //n = (nocens * n) + (cens * Nvec);
  
  for(int i = 0;i <  N.size(); ++i){
    cens[i] = 0;
    if (n[i] > Nvec[i]) cens[i] = 1;
  }
  
  nocens = 1 - cens;
  
  arma::vec outvec1(N.size());
  arma::vec outvec2(N.size());
  for(int i = 0;i <  N.size(); ++i){
    if(minnumstrat < Nvec[i]){outvec[i]= minnumstrat;}else{outvec[i]=Nvec[i];};
    outvec1[i]=(nocens[i] * n[i]);
  }
  //cout << "outvec1: " << outvec1 << "\n  ";
  for(int i = 0;i <  N.size(); ++i){
    outvec2[i]=(cens[i] * Nvec[i]);
  }
  //cout << "outvec2: " << outvec2 << "\n  ";
  n = outvec1 + outvec2;
  return n;
  
  
}
