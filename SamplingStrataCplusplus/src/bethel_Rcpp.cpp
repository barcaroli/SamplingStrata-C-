//-------------------------------------------------------
// C++ script implementing Bethel algorithm
// Author: Giulio Barcaroli
// June 2021
//-------------------------------------------------------

#include <Rcpp.h>
#include <iostream>

using std::cout;
using std::cin;
using std::endl;

using namespace Rcpp;


//--------------  Routine to select variables
// [[Rcpp::export]]
NumericMatrix select_variables(DataFrame dati, 
                           std::string prefix, 
                           int nvar)    {
  StringVector cols (nvar);
  for ( int i=0; i<nvar; i++ ) {
    std::ostringstream str1; 
    str1 << i+1; 
    std::string str2 = str1.str();
    std::string s1(prefix);
    std::string y = s1 + str2;
    cols[i] = y;
    //    std::cout << "cols: " << cols << std::endl;
  }
  DataFrame subset = dati[cols];
  NumericMatrix mat = internal::convert_using_rfunction(subset, "as.matrix");  
  return(mat);
}


//--------------  Routine to create disjoint matrix
// [[Rcpp::export]]
IntegerMatrix disjoint(NumericVector dom){ 
  double nvalues = dom[which_max(dom)];
  int nstrat = dom.length();
  IntegerMatrix disj(nstrat,nvalues);
  for ( int i = 0; i < nvalues; i++ ) {
    for ( int j = 0; j < nstrat; j++ ) {
      //cout << "dom" << dom[j] << endl;
      //cout << "i" << i << endl;
      if (dom[j] == i+1) disj(j,i) = 1;
      if (dom[j] != i+1) disj(j,i) = 0;
    }
  }
  return(disj);
}

//--------------  Routine for matrices product
// [[Rcpp::export]]
NumericMatrix m_s(IntegerMatrix disj,
                  NumericMatrix mat) {
  int nc = disj.ncol();
  //cout << "nc " << nc << endl;
  int nstrat = disj.nrow();
  //cout << "nstrat " << nstrat << endl;
  int nvar = mat.ncol();
  //cout << "nvar " << nvar << endl;
  NumericMatrix m(nstrat,nvar);
  NumericMatrix m1(nstrat,nvar);
  NumericMatrix m2(nstrat,nvar);
  for (int k = 0; k < nc; k++) {
    //cout << "k " << k << endl;
    for (int i = 0; i < nstrat; i++)  {
      for (int j = 0; j < nvar; j++) {
        m1(i,j) = mat(i,j) * disj(i,k) ;
      }
    }
    m = cbind(m,m1); 
  }
  //cout << "m " << m << endl;
  int t = m.ncol();
  NumericMatrix::Sub    sub = m( Range(0,nstrat-1) , Range(nvar,t-1) );
  return(sub);
}

//--------------  rowSums
// [[Rcpp::export]]
NumericVector rowSums_Rcpp(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nr);
  for (int i = 0; i < nr; i++) {
    double sum = 0.0;
    for (int j = 0; j < nc; j++) {
      sum += x(i, j);
    }
    ans[i] = sum;
  }
  return ans;
}

//--------------  colSums
// [[Rcpp::export]]
NumericVector colSums_Rcpp(const NumericMatrix& x) {
  int nr = x.nrow(), nc = x.ncol();
  NumericVector ans(nc);
  for (int j = 0; j < nc; j++) {
    double sum = 0.0;
    for (int i = 0; i < nr; i++) {
      sum += x(i, j);
    }
    ans[j] = sum;
  }
  return ans;
}

//------------------------ Routine linearize cv's
// [[Rcpp::export]]
NumericMatrix cv_Rcpp(DataFrame errors, 
                      int ndom, 
                      int nvar)    {
  NumericMatrix cv1(1,nvar);
  NumericMatrix cvx(1,nvar);
  for ( int k = 0; k < ndom; k++ ) {
    cvx = select_variables(errors,"CV",nvar);
    cv1 = cbind(cvx,cv1);
  }
  NumericMatrix::Sub    cv = cv1( Range(0,0) , Range(0,(nvar*ndom-1)) );
  return(cv);
}

//------------------------ Routine crea_a
// [[Rcpp::export]]
NumericMatrix crea_a(NumericMatrix& m,
                     NumericMatrix& s,
                     NumericVector& nocens,
                     NumericVector& N,
                     NumericVector& cv,
                     double& epsilon){
   int nr = m.nrow(), nc = m.ncol();
  NumericMatrix numA(nr,nc);
  for ( int i = 0; i < s.nrow(); i++ ) {
    for ( int j = 0; j < s.ncol(); j++ ) {
      numA(i,j) = N(i) * N(i) * s(i,j)*s(i,j) * nocens(i);
    }
  }
  
  NumericMatrix y(nr,nc);
  NumericVector denA1(nc); 

  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      y(i,j) = N(i)*m(i,j);
    }
  }

  y = transpose(y);

  for (int j = 0; j < nr; j++) {
    for (int i = 0; i < nc; i++) {
      y(i,j) = y(i,j) * cv(i);
    }
  }

  y = transpose(y);

  denA1 = colSums_Rcpp(y);
  denA1 = pow(denA1,2);
  
  NumericMatrix w(s.nrow(),s.ncol());
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      w(i,j) = N(i)*pow(s(i,j),2)*nocens(i);
    }
  }
  NumericVector denA2 = colSums(w);
  
  NumericVector denA = denA1 + denA2 + epsilon;
  
  NumericMatrix z(nc,nr);
  z = transpose(numA);
  
  for (int j = 0; j < nr; j++) {
    for (int i = 0; i < nc; i++) {
      z(i,j) = z(i,j) / denA(i);
    }
  }
  
  NumericMatrix a(nr,nc);
  a = transpose(z);

  return(a);
}




//--------------  Chromy algorithm
// [[Rcpp::export]]
NumericVector chromy_Rcpp(NumericMatrix a,
                          double alfatot, 
                          double diff, 
                          int iter, 
                          NumericVector alfa, 
                          NumericVector alfanext, 
                          NumericVector x,
                          NumericVector cost,
                          int nvar,
                          bool realAllocation
) {
  int maxiter = 200; 
  
  
  double epsilon = 1e-11;
  int nr=a.nrow();
  int nc=nvar;
  //std::cout << "nr \n" << nr << std::endl;
  //std::cout << "nc \n" << nc << std::endl;
  
  NumericVector n(nr);
  
  while (diff > epsilon && iter < maxiter) {
    
    
    //den1 <- sqrt(rowSums(t(t(a) * c(alfa))))
    NumericMatrix b(nc,nr);
    b = transpose(a);
    
    for (int j = 0; j < nr; j++) {
      for (int i = 0; i < nc; i++) {
        b(i, j) = b(i, j) * alfa(i);
      }
    }
    
    NumericMatrix c(nr,nc);
    c = transpose(b);
    
    NumericVector den1(nr);
    den1 = pow(rowSums_Rcpp(c),0.5);
    //std::cout << "den1 \n" << den1 << std::endl;
    
    //den2 <- sum(sqrt(rowSums(t(t(a * cost) * c(alfa)))))
    
    for (int j = 0; j < nc; j++) {
      for (int i = 0; i < nr; i++) {
        c(i, j) = a(i, j) * cost(i);
      }
    }
    b = transpose(c);
    
    for (int j = 0; j < nr; j++) {
      for (int i = 0; i < nc; i++) {
        b(i, j) = b(i, j) * alfa(i);
      }
    }
    
    c = transpose(b);
    
    NumericVector d(nr);
    d = pow(rowSums_Rcpp(c),0.5);
    
    double den2;
    den2 = 0;
    
    for (int j = 0; j < nr; j++) {
      den2 = den2 + d(j);
    }  
    //std::cout << "den2 \n" << den2 << std::endl;
    
    //x <- sqrt(cost)/(den1 * den2 + epsilon)
    
    for (int j = 0; j < nr; j++) {
      x(j) = pow(cost(j),0.5) / (den1(j)*den2+epsilon);
    } 
    //std::cout << "x \n" << x << std::endl;
    
    //alfatot <- sum(c(alfa) * (t(a) %*% x)^2)
    
    b = transpose(a);
    Function f1("%*%"); 
    NumericVector e(nc);
    e = f1(b,x);
    for (int j = 0; j < nc; j++) {
      e(j) = alfa(j) * pow(e(j),2);
    }
    alfatot = 0;
    for (int j = 0; j < nc; j++) {
      alfatot += e(j);
    } 
    
    //alfatot[alfatot == 0] <- epsilon
    
    if (alfatot == 0) alfatot = epsilon;
    //cout<<"alfatot \n"<< alfatot <<endl;
    
    //alfanext <- c(alfa) * (t(a) %*% x)^2/alfatot
    
    alfanext = e/alfatot;
    //cout<<"alfanext \n"<< alfanext <<endl;
    //cout<<"alfa \n"<< alfa <<endl;
    //std::cout << "alfanext \n" << alfanext << std::endl;
    
    //diff <- max(abs(alfanext - alfa))
    
    
    NumericVector alfadiff(nc);
    for (int j = 0; j < nc; j++) {
      alfadiff(j) = alfanext(j) - alfa(j);
    }
    
    //std::cout << "alfadiff \n" << alfadiff << std::endl;
    
    double max;
    max=0;
    int k;
    k=1;
    
    for (int j = 0; j < nc; j++) {
      if (alfadiff(j) < 0) alfadiff(j) = alfadiff(j) * -1;
      if (alfadiff(j) > max) {
        max = alfadiff(j);
        k = j;
      }
    }  
    //cout<<"alfadiff \n"<< alfadiff <<endl;
    diff = alfadiff(k);
    //cout<<"diff \n"<< diff <<endl;
    
    for (int j = 0; j < nc; j++) {
      alfa(j) = alfanext(j);
    }
    
    //if (realAllocation == 1) n = 1/x;
    //if (realAllocation == 0) n = ceil(1/x);
    //cout<<"n \n"<< n <<endl;
  }
  
  if (realAllocation == 1) n = 1/x;
  if (realAllocation == 0) n = ceil(1/x);
  
  return(n);
}

//-------------------------- check_n
// [[Rcpp::export]]
NumericVector check_n(NumericVector n,
                      NumericVector N,
                      int minnumstrat) {
  int nstrat = n.length();
  NumericVector n1(nstrat);
  
  for (int i = 0; i < nstrat; i++) {
    n1(i) = n(i);
    if (n(i) > N(i)) n1(i) = N(i);
    if (n(i) < minnumstrat) {
      if (N(i) >= minnumstrat) n1(i) = minnumstrat;
      if (N(i) < minnumstrat) n1(i) = N(i);
    }
  }
  return(n1);
}

//--------------  Bethel allocation
// [[Rcpp::export]]
NumericVector bethel_Rcpp(DataFrame strata,
                          DataFrame errors,
                          int minnumstrat,
                          bool realAllocation = 1) {

  //int maxiter = 200; 
  //int maxiter1 = 25;
  double epsilon = 1e-11;

  int nstrat = strata.nrows();
  //cout << "nstrat: " << nstrat << endl;
  StringVector errors_names = errors.names();
  //cout << "cv names: " << errors_names << endl;
  int nvar = errors_names.length() - 2;
  bool isPresent = (std::find(errors_names.begin(), errors_names.end(), "domainvalue") != errors_names.end());
  if (isPresent == FALSE) nvar = errors_names.length() - 1;
  //cout << "cv: " << nvar << endl;
  

  NumericMatrix med = select_variables(strata,"M",nvar);

  NumericMatrix esse = select_variables(strata,"S",nvar);

  NumericVector N = strata["N"];
  NumericVector cens = strata["CENS"];
  /*
  for ( int i = 0; i < nstrat; i++ ) {
    if (cens[i] < minnumstrat) cens[i] = 1;
  }
  */
  NumericVector nocens = 1 - cens;
  NumericVector cost = strata["COST"];
  
  std::string nomdom = "DOM1";
  NumericVector dom = strata["DOM1"];
  int ndom = dom[which_max(dom)];
  //cout << "ndom: " << ndom << endl;
  
  IntegerMatrix disj = disjoint(dom);
  
  NumericMatrix m = m_s(disj,med);
  //cout << "m: " << m << endl;
  NumericMatrix s = m_s(disj,esse);
  //cout << "s: " << s << endl;
  
  NumericMatrix cv(1,nvar);
  
  cv = cv_Rcpp(errors,ndom,nvar);
  //cout << "cv: " << cv << endl;
  
  nvar = cv.ncol();
  //cout << "nvar: " << nvar << endl;

  NumericMatrix a(nstrat,nvar);
  //std::cout << "crea " << a << std::endl;
  a = crea_a(m,s,nocens,N,cv,epsilon);
  //cout << "a \n " << a << endl;
  
  double alfatot;
  alfatot = 0;
  double diff;
  diff= 999;
  int iter;
  iter = 0;
  /*
  bool realAllocation;
  realAllocation = 0;
  */
  
  NumericVector alfa(nvar);
  NumericVector alfanext(nvar);

  for (int j = 0; j < nvar; j++) {
    alfa(j) = 1.0/nvar;
  }

  NumericVector x(nstrat);  
  for (int i = 0; i < nstrat; i++) {
    x(i) = 0.1;
    //cout << "x(i): " << x(i) << endl;
  }

  /*
  cout << "alfatot: " << alfatot << endl;
  cout << "diff: " << diff << endl;
  cout << "iter: " << iter << endl;
  cout << "alfa: " << alfa << endl;
  cout << "alfanext: " << alfanext << endl;
  cout << "x: " << x << endl;
  cout << "cost: " << cost << endl;
  cout << "nvar: " << nvar << endl;
  cout << "realAllocation: " << realAllocation << endl;
   */
  
  NumericVector n(nstrat);
  n = chromy_Rcpp(a,alfatot,diff,iter,alfa,alfanext,x,cost,nvar,realAllocation);
  //cout << "n: " << n << endl;
  
  int contx;
  contx = 0;
  for (int i = 0; i < nstrat; i++) {
    if (n(i) > N(i)) {
      contx += 1;
      cens(i) = 1;
      nocens(i) = 0;
    }
  }
  
  //NumericVector n2(nstrat);
  n = check_n(n,N,minnumstrat);
  //cout << "n: " << n << endl;
  
  int iter1 = 0;
  int maxiter1 = 25;
  
  while (contx > 0 && iter1 < maxiter1) {
    iter1 = iter1 + 1;
    a = crea_a(m,s,nocens,N,cv,epsilon);
    n = chromy_Rcpp(a,alfatot,diff,iter,alfa,alfanext,x,cost,nvar,realAllocation);
    contx = 0;
    for (int i = 0; i < nstrat; i++) {
      if (n(i) > N(i)) {
        contx += 1;
        cens(i) = 1;
        nocens(i) = 0;
        n(i) = N(i);
      }
    }
    n = check_n(n,N,minnumstrat);
    //cout << "n: " << n << endl;
  }
  
  //n <- (nocens * n) + (cens * N)

  for (int i = 0; i < nstrat; i++) {
       n(i) = n(i)*nocens(i) + N(i)*cens(i);
    }

  
   return(n);
}



