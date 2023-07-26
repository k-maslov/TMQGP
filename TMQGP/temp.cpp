#include "temp.h"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <complex>

#include <omp.h>
#include <math.h>
#include <vector>




double full_int(double omega, double m, double eps){
    Int_gsl_adaptive inte;
    double resu, erro;
    gsl_set_error_handler_off();
    funct h_lambda = [&] (double p){
        funct h_lambda1 = [&] (double p0){
         //   double sign_p0 = (p0 >= 0) ? 1.0 : -1.0;



//            return  imag(-4*3*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(4*3.14*3.14*3.14*(p0*p0 - p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)*(p0*p0 + p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)));
  //          };        
//           return imag(-4*3*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(4*3.14*3.14*3.14*(( (p0-omega/2.0)*(p0-omega/2.0) - p*p - m*m) + complex<double> (0,1)*eps)*((p0+(omega/2))*(p0+(omega/2)) - p*p - m*m + complex<double> (0,1)*eps)));
//        };
//        IntGSL<funct> inte2;



        return  imag(-4*3*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(4*3.14*3.14*3.14*(p0*p0 - p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)*(p0*p0 + p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)));
        };


        Int_gsl_adaptive inte2;
        double resu1, erro1;
        inte2.integrate(&h_lambda1, -1.641, 1.641, resu1, erro1); // limit = 1.641 for m= 0.5
        return resu1;
    };
    inte.integrate(&h_lambda, 0, 0.695, resu, erro);
    return resu;
}

double p0_integral(double omega, double p, double m, double eps){
    funct i_lambda_1 = [&] (double p0){
        return imag(-complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(((p0-(omega/2.0))*(p0-(omega/2.0)) - p*p - m*m + complex<double> (0,1)*eps)*((p0+(omega/2))*(p0+(omega/2)) - p*p - m*m + complex<double> (0,1)*eps)));
    };

    Int_gsl_adaptive integ1;
    gsl_set_error_handler_off();
    double res1, err1;
    integ1.integrate(&i_lambda_1, -1.641, 1.641, res1, err1);
    return res1;
}

// iNTEGRAND FOR VACUUM POLARIZATIOION
double Integrand(double omega, double p0, double p, double m, double eps){
   return imag(-12*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(4*3.14*3.14*3.14*(p0*p0 - p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)*(p0*p0 + p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)));
}


double full_int_using_interpo(double omega, double m, Interpolator2D & ImG){
    Int_gsl_adaptive inte_p0;
    double resu_p0, erro_p0;
    //gsl_set_error_handler_off();
    funct p0_lambda = [&] (double p){
        funct p_lambda = [&] (double p0){
            return imag(-12 * ImG(p, p0 + omega/2) * ImG(p, p0 - omega/2)* p*p * complex<double> (0,1)); //* (p*p - p0*p0 +(omega*omega/4) +m*m)*(tanh((p0-omega)/2)/2*T + tanh((p0+omega)/2*T))/ (4*3.14*3.14*3.14));
        };
        Int_gsl_adaptive inte_p;
        double resu_p, erro_p;
        inte_p.integrate(&p_lambda, 0, 0.695, resu_p, erro_p); // limit = 1.641 for m= 0.5
        return resu_p;

    };
    inte_p0.integrate(&p0_lambda, -1.641, 1.641, resu_p0, erro_p0);
    return resu_p0;

}

double p_integral(double omega, double p0, double m, double eps){
    funct i_lambda_p = [&] (double p){
        return imag(-p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/((( (p0-omega/2.0)*(p0-omega/2.0) - p*p - m*m) + complex<double> (0,1)*eps)*((p0+(omega/2))*(p0+(omega/2)) - p*p - m*m + complex<double> (0,1)*eps)));
    };

    Int_gsl_adaptive integ_p;
    gsl_set_error_handler_off();
    double res_p, err_p;
    integ_p.integrate(&i_lambda_p, 0, 0.695, res_p, err_p);
    return res_p;
}


double FInt_no_inter(double omega, double p, double m, double eps){

    gsl_set_error_handler_off();

    Int_gsl_adaptive integ_a;
    Int_gsl_adaptive integ_b;
 

    funct h_lambda_p1 = [&] (double p_1){
        funct h_lambda_x = [&] (double omega_1){
          
            double sign_o = (omega - omega_1 >= 0) ? 1.0 : -1.0;
            double sign_o1 = (omega_1 >= 0) ? 1.0 : -1.0;
         
            return imag(p_1 * p_1 /(4*3.14*3.14*3.14*((omega - omega_1)*(omega -omega_1) - (p - p_1)*(p - p_1)- m*m + complex<double> (0,1)*eps*sign_o) * ((omega_1 * omega_1) - (p_1 * p_1) - m*m + complex<double> (0,1)*eps*sign_o1)));
            //return p_1 * p_1 /(((omega -omega_1)*(omega -omega_1) - (p - p_1)*(p - p_1)- m*m + complex<double> (0,1)*eps) * ((omega_1 * omega_1) - (p_1 * p_1) - m*m + complex<double> (0,1)*eps))); 
        };
        double res_a, err_a;
        integ_b.integrate(&h_lambda_x,-1.641, 1.641, res_a, err_a);
        return res_a;

    };
    gsl_set_error_handler_off();
    
    double res_c, err_c;
    integ_a.integrate(&h_lambda_p1, 0, 0.695, res_c, err_c);
    return res_c;
}


// // // ////
////////
/////////
 // / / // /// // 
 // THERMAL CALCULATIONS FUNCTIONS:

double FInt_t(double omega, double p, double T, Interpolator2D & ImG){

    gsl_set_error_handler_off();

    // IntGSL<std::function<double(double)>> integ;
    Int_gsl_adaptive integ_t;
    Int_gsl_adaptive integ_1_t;
    //IntGSL<funct> integ_1;
    double res, err;
    
    funct i_lambda_p1_t = [&] (double p_1){
        funct i_lambda_t = [&] (double omega_1){
            return ImG(p - p_1, omega - omega_1)* ImG(p_1, omega_1) * (tanh((omega-omega_1)/(2*T)) + tanh(omega_1/(2*T)));
            };
    integ_t.integrate(&i_lambda_t,-1.641, 1.641, res, err);
    return res;
        };

    gsl_set_error_handler_off();
    
    double res_p1, err_p1;
    integ_1_t.integrate(&i_lambda_p1_t, 0, 0.695, res_p1, err_p1);
    return res_p1;
    }
//
//
// Plotting the integral 6.31 from Hatsuda 


double F_p(double omega,double q,double T,double m,double eps, double lambda){
    Int_gsl_adaptive int_t;
    double resu_t, erro_t;
    gsl_set_error_handler_off();
    funct h_lambda_t = [&] (double p){
        funct h_lambda1_t = [&] (double p0){
            int N_c = 3;
            double sign_p01 = (p0 + omega >= 0) ? 1.0 : -1.0;
            double sign_p02 = (p0 - omega >= 0) ? 1.0 : -1.0;
            double sign_p03 = (p0  >= 0) ? 1.0 : -1.0;
            complex<double> coeff1 = 4 * 4 * 3 * p * p * tanh(p0/2/T) *  (m*m - p0*p0 + p*p - p0*omega)/8/M_PI/M_PI/M_PI ;
            complex<double> coeff2 = 4 * 4 * 3 * p * p * tanh(p0/2/T) *  (m*m - p0*p0 + p*p + p0*omega)/8/M_PI/M_PI/M_PI ; //10.08 * 10.08 *
            complex<double> green1 = complex<double>(1,0)/((p0 + omega)*(p0 + omega) - (p+q)*(p+q) - m*m + complex<double>(0,1) *eps* sign_p01) ; 
            complex<double> green2 = complex<double>(1,0)/(p0*p0 - p*p - m*m + complex<double>(0,1) * eps * sign_p03);
            complex<double> green3 = complex<double>(1,0)/((p0 - omega)*(p0 - omega) - (p - q)*(p - q) - m*m - complex<double>(0,1) *eps * sign_p02);

           // return real(coeff1 * green1 * imag(green2) + coeff2*imag(green2) * green3);
            return imag(coeff1 * green1 * imag(green2) + coeff2*imag(green2) * green3);
        };
        Int_gsl_adaptive inte2_t;
        
        double result, error;
        inte2_t.integrate(&h_lambda1_t, -sqrt(4 * (lambda*lambda + m*m)),sqrt(4 * (lambda*lambda + m*m)) , result, error); // limit = 1.641 for m= 0.5
        return result;
    };
    int_t.integrate(&h_lambda_t, 0, lambda, resu_t, erro_t);
    return resu_t;
}

double F_s(double omega,double q,double T,double m,double eps, double lambda){
    Int_gsl_adaptive int_t;
    double resu_t, erro_t;
    gsl_set_error_handler_off();
    funct h_lambda_t = [&] (double p){
        funct h_lambda1_t = [&] (double p0){
            double sign_p01 = (p0 + omega >= 0) ? 1.0 : -1.0;
            double sign_p02 = (p0 - omega >= 0) ? 1.0 : -1.0;
            double sign_p03 = (p0  >= 0) ? 1.0 : -1.0;
            complex<double> coeff1 = 4 * 4 * 3 * p * p * tanh(p0/2/T) *  (m*m + p0*p0 - p*p + p0*omega)/8/M_PI/M_PI/M_PI ;
            complex<double> coeff2 = 4 * 4 * 3 * p * p * tanh(p0/2/T) *  (m*m + p0*p0 - p*p - p0*omega)/8/M_PI/M_PI/M_PI ; //10.08 * 10.08 *
            complex<double> green1 = complex<double>(1,0)/((p0 + omega)*(p0 + omega) - (p)*(p) - m*m + complex<double>(0,1) *eps* sign_p01) ; 
            complex<double> green2 = complex<double>(1,0)/(p0*p0 - p*p - m*m + complex<double>(0,1) * eps * sign_p03);
            complex<double> green3 = complex<double>(1,0)/((p0 - omega)*(p0 - omega) - (p)*(p) - m*m - complex<double>(0,1) *eps * sign_p02);

            //return real(coeff1*(real(green1) * imag(green2) + imag(green2) * real(green3))); 
            return imag(coeff1 * green1 * imag(green2) + coeff2*imag(green2) * green3);
        };
        Int_gsl_adaptive inte2_t;
        double result, error;
        inte2_t.integrate(&h_lambda1_t, -sqrt(4 * (lambda*lambda + m*m)),sqrt(4 * (lambda*lambda + m*m)) , result, error); // limit = 1.641 for m= 0.5
        return result;
    };
    int_t.integrate(&h_lambda_t, 0, lambda, resu_t, erro_t);
    return resu_t;
}
// double Fp_real(double omega, double q, double T, double m , double eps){
//     gsl_set_error_handler_off();
//     double res, err;

//     funct i_func = [&] (double omega1) -> double{
        
//     };

//     inter cauchy.integrate(&i_func, -3, 3, omega, res, err );
//     return res;
// }
double Re_meson(double omega, double lambda4,Interpolator & iImS){
   // gsl_set_error_handler_off();
    double res, err;

    funct i_func = [&] (double omega1) -> double{
        return iImS(omega1) / M_PI;
    };
    Int_gsl_cauchy inte_c;
    inte_c.integrate(&i_func, -lambda4 , lambda4, omega, res, err);
    
    return res;
}


double imag_pi_inter1( double omega,double T,double m, double lambda, Interpolator2D & ImG){
    //Int_gsl_adaptive int_t;
    double resu_t, erro_t;
    //gsl_set_error_handler_off();
    funct h_lambda_t = [&] (double p0){
        funct h_lambda1_t = [&] (double p){

            double coeff1 = 4 * 4 * 3 * p * p * tanh(p0/2/T) *  (m*m - p0*p0 + p*p - p0*omega)/8/M_PI/M_PI/M_PI ;
            double coeff2 = 4 * 4 * 3 * p * p * tanh(p0/2/T) *  (m*m - p0*p0 + p*p + p0*omega)/8/M_PI/M_PI/M_PI ;
            
            return coeff1 * ImG(p, p0 + omega) * (ImG(p, p0))  + coeff2 * ImG(p , p0) * ImG(p, omega - p0);
        };
        //Int_gsl_adaptive inte2_t;
        double result, error;
       // gsl_set_error_handler_off();
        integ_kt.integrate(&h_lambda1_t, 0 , lambda , result, error); // limit = 1.641 for m= 0.5
        return result;
    };
    integ_Et.integrate(&h_lambda_t, -sqrt(4 * (lambda*lambda + m*m)), sqrt(4 * (lambda*lambda + m*m)), resu_t, erro_t);
    return resu_t;  
}



// double Sigma_re_plus1(double p0, double T, double m, double lambda, Interpolator2D & ImG){
//     double res, err;
//     Int_gsl_adaptive int_t;
//     funct h_lambda_t = [&] (double q){
    
//         double Eq = sqrt(q*q + m*m);
//         double coeffp = (1/tanh((-Eq + p0)/2/T) + tanh(Eq/2/T))*(1 + m/Eq);

//         return q*q*4*M_PI * (M_PI/2) * (1/pow(2*M_PI, 4)) * ImG(q , -Eq + p0)*coeffp; 
//         };

            
   
// integ_kt.integrate(&h_lambda_t,-0.631,0.631 /*-sqrt(4 * (lambda*lambda + m*m)), sqrt(4 * (lambda*lambda + m*m))*/, res, err);
// return res;  
// }



double Sigma_re_plus1(double p0, double T, double mf, double mb, Interpolator2D & ImG){
    double res, err;
    gsl_set_error_handler_off();
    Int_gsl_adaptive int_t;
    funct h_lambda_t = [&] (double p){ 

        double Ep = sqrt(p*p + mf*mf);

        double sign_p01 = (p0 + Ep >= 0) ? 1.0 : -1.0;
        double sign_p02 = (p0 - Ep >= 0) ? 1.0 : -1.0;

        double coeffp = p*p*(1 - mf/Ep)/8/M_PI/M_PI;
        coeffp *= (1/tanh((p0 - Ep)/2/T) + tanh(Ep/2/T)) ;
        coeffp *= ImG(p, p0 - Ep);
        //coeffp *= imag(complex<double>(1,0)/((p0-Ep)*(p0-Ep) - p*p - mb*mb + complex<double>(0,1)*eps*sign_p02));

        return coeffp;
        };
 
int_t.integrate(&h_lambda_t, 0 , 2  /*-sqrt(4 * (lambda*lambda + m*m)), sqrt(4 * (lambda*lambda + m*m))*/, res, err);
return res;  
}

double Sigma_re_minus1(double p0, double T, double mf, double mb, Interpolator2D & ImG){
    double res, err;
    gsl_set_error_handler_off();
    Int_gsl_adaptive int_t;
    funct h_lambda_t = [&] (double p){ 
    
        double Ep = sqrt(p*p + mf*mf);

        double sign_p01 = (p0 + Ep >= 0) ? 1.0 : -1.0;
        double sign_p02 = (p0 - Ep >= 0) ? 1.0 : -1.0;

        double coeffp = p*p*(1 + mf/Ep)/8/M_PI/M_PI;
        coeffp *= (1/tanh((p0 + Ep)/2/T) + tanh(-Ep/2/T)) ;
        coeffp *= ImG(p , p0 + Ep);
        //coeffp *= imag(complex<double>(1,0)/((p0+Ep)*(p0+Ep) - p*p - mb*mb + complex<double>(0,1)*eps*sign_p01));

        return coeffp;
        };
int_t.integrate(&h_lambda_t, 0 , 2, res, err);
return res;   
}

// QUARK LOOPS

/* This is for q = 0*/
complex<double> G0(double p0, double p, double m, double eps){
    double sign_p0 = (p0 >= 0) ? 1.0 : -1.0;
    return (p0/(p0*p0 - p*p - m*m + complex<double>(0,1)* eps * sign_p0));
}

complex<double> Gv(double p0, double p, double m, double eps){
    double sign_p0 = (p0 >= 0) ? 1.0 : -1.0;
    return (complex<double>(p,0)/(p0*p0 - p*p - m*m + complex<double>(0,1)* eps * sign_p0));
}

complex<double> Gm(double p0, double p, double m, double eps){
    double sign_p0 = (p0 >= 0) ? 1.0 : -1.0;
    return (m/(p0*p0 - p*p - m*m + complex<double>(0,1)* eps * sign_p0));
}


double Impi(double p, double omega, double T, double m, double eps){ /*q = 0*/
    gsl_set_error_handler_off();

    double err, res, res1, err1;

    funct lambda2 = [&] (double z){
        funct lambda3 = [&] (double x) {

            double coeff = p * p * imag(G0(z, p, m, eps)) *  imag(G0(z - omega, p, m, eps));
        
            coeff += p * p * x * imag(Gv(z, p, m, eps)) *  imag(Gv(z - omega, p, m, eps)); /* This term is zero for q = 0*/

            coeff += - p * p * imag(Gm(z, p, m, eps)) *  imag(Gm(z - omega, p, m, eps));

            coeff *= (tanh((z - omega)/2/T) - tanh(z/2/T))/2/M_PI/M_PI/M_PI;

            return coeff;

    };

        integ_Et.integrate(&lambda3 , -1, 1, res, err);
        return res;
    };
    integ_kt.integrate(&lambda2 , -1.651, 1.651, res1, err1);
    return res1;
}

double Impi_int(double omega, double T, double m, double eps){
    gsl_set_error_handler_off();

    double res,err;
    Int_gsl_adaptive int_t;
    funct lambda = [&] (double p){

        return Impi(p, omega, T, m, eps);
    };
    int_t.integrate(&lambda, 0, 1.651, res, err);
    return res;

}

// HARADA NEMOTO


double Impi2(double p, double q, double omega, double T, double m, double eps){
    gsl_set_error_handler_off();

    double err, res, res1, err1;

    funct lambda2 = [&] (double z){
        funct lambda3 = [&] (double k) {

            double coeff = p * k * imag(G0(z, p, m, eps)) *  imag(G0(z - omega, k, m, eps));

            coeff += -(q * q - p * p - k * k ) *  imag(Gv(z, p, m, eps)) * imag(Gv(z - omega, k, m , eps))/2;
            //coeff += -(q * q + p * p - k * k ) * (1 - p/q) * imag(Gv(z, p, m, eps))/2 * imag(Gv(z - omega, k, m , eps));

            coeff += - p * k * imag(Gm(z, p, m, eps)) *  imag(Gm(z - omega, k, m, eps));

            coeff *= p *(tanh((z - omega)/2/T) - tanh(z/2/T))/2/M_PI/M_PI/M_PI/(q);
            //coeff *= p *(tanh((z - omega)/2/T) - tanh(z/2/T))/2/M_PI/M_PI/M_PI/(p - q);

            return coeff;

    };

        integ_Et.integrate(&lambda3 , abs(p - q), p + q, res, err);
        return res;
    };
    integ_kt.integrate(&lambda2 , -2, 2, res1, err1);
    return res1;
}

double Impi_int2(double q , double omega, double T, double m, double eps){
    gsl_set_error_handler_off();

    double res,err;
    Int_gsl_adaptive int_t;
    funct lambda = [&] (double p){

        return Impi2(p, q, omega, T, m, eps);
    };
    integ_kt.integrate(&lambda, 0, 1, res, err);
    return res;

}