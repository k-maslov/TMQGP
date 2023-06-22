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
    gsl_set_error_handler_off();
    funct p0_lambda = [&] (double p){
        funct p_lambda = [&] (double p0){
            return imag(-12 * ImG(p, p0 + omega/2) * ImG(p, p0 - omega/2)* p*p * complex<double> (0,1)); //* (p*p - p0*p0 +(omega*omega/4) +m*m)*(tanh((p0-omega)/2)/2*T + tanh((p0+omega)/2*T))/ (4*3.14*3.14*3.14));
        };
        Int_gsl_adaptive inte_p;
        double resu_p, erro_p;
        inte_p.integrate(&p_lambda, -1.641, 1.641, resu_p, erro_p); // limit = 1.641 for m= 0.5
        return resu_p;

    };
    inte_p0.integrate(&p0_lambda, 0, 0.695, resu_p0, erro_p0);
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


double Fint_integrand(double omega,double omega_1, double p, double p_1, double T, double m, double eps){
    double sign_o = (omega - omega_1 >= 0) ? 1.0 : -1.0;
    double sign_o1 = (omega_1 >= 0) ? 1.0 : -1.0;
         
    return imag((tanh((omega-omega_1)/(2*T)) + tanh(omega_1/(2*T)))/(16*3.14*3.14*3.14*3.14*((omega -omega_1)*(omega -omega_1) - (p - p_1)*(p - p_1)- m*m + complex<double> (0,1)*eps*sign_o) * ((omega_1 * omega_1) - (p_1 * p_1) - m*m + complex<double> (0,1)*eps*sign_o1)));
}


double FInt_no_inter_t(double omega, double p, double T, double m, double eps){

    gsl_set_error_handler_off();

    Int_gsl_adaptive integ_a_t;
    Int_gsl_adaptive integ_b_t;
 

    funct h_lambda_p1_t = [&] (double p_1){
        funct h_lambda_x_t = [&] (double omega_1){
          
            double sign_o = (omega - omega_1 >= 0) ? 1.0 : -1.0;
            double sign_o1 = (omega_1 >= 0) ? 1.0 : -1.0;
         
            return imag(p_1 * p_1 *(tanh((omega-omega_1)/(2*T)) + tanh(omega_1/(2*T)))/(4*3.14*3.14*3.14*((omega -omega_1)*(omega -omega_1) - (p - p_1)*(p - p_1)- m*m + complex<double> (0,1)*eps*sign_o) * ((omega_1 * omega_1) - (p_1 * p_1) - m*m + complex<double> (0,1)*eps*sign_o1)));
            //return (tanh((omega-omega_1)/(2*T)) + tanh(omega_1/(2*T)))/(((omega -omega_1)*(omega -omega_1) - (p - p_1)*(p - p_1)- m*m + complex<double> (0,1)*eps) * ((omega_1 * omega_1) - (p_1 * p_1) - m*m + complex<double> (0,1)*eps))); 
        };
        double res_a, err_a;
        integ_b_t.integrate(&h_lambda_x_t,-1.641, 1.641, res_a, err_a);
        return res_a;

    };
    gsl_set_error_handler_off();
    
    double res_c, err_c;
    integ_a_t.integrate(&h_lambda_p1_t, 0, 0.695, res_c, err_c);
    return res_c;
}

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


double full_int_t(double omega, double T, double m, double eps){
    Int_gsl_adaptive inte_t;
    double resu_t, erro_t;
    gsl_set_error_handler_off();
    funct h_lambda_t = [&] (double p){
        funct h_lambda1_t = [&] (double p0){
            double sign_p0 = (p0 + omega >= 0) ? 1.0 : -1.0;
            double sign_p0_omega = (p0 - omega >= 0) ? 1.0 : -1.0;
//            return  imag(-4*3*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(4*3.14*3.14*3.14*(p0*p0 - p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)*(p0*p0 + p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps)));
  //          };        
//           return imag(-4*3*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)/(4*3.14*3.14*3.14*(( (p0-omega/2.0)*(p0-omega/2.0) - p*p - m*m) + complex<double> (0,1)*eps)*((p0+(omega/2))*(p0+(omega/2)) - p*p - m*m + complex<double> (0,1)*eps)));
//        };
//        IntGSL<funct> inte2;
        return  imag(-4*2*3*p*p* complex<double> (0,1)*(p*p - p0*p0 +(omega*omega/4) +m*m)*(tanh((p0-(omega/2))/2/T) + tanh((p0+omega/2)/2/T))/(4*M_PI*M_PI*M_PI*(p0*p0 - p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps*sign_p0_omega)*(p0*p0 + p0*omega + omega*omega/4 - p*p - m*m + complex<double> (0,1)*eps*sign_p0)));
        };
        Int_gsl_adaptive inte2_t;
        double resu1_t, erro1_t;
        inte2_t.integrate(&h_lambda1_t, -1.641, 1.641, resu1_t, erro1_t); // limit = 1.641 for m= 0.5
        return resu1_t;
    };
    inte_t.integrate(&h_lambda_t, 0, 0.695, resu_t, erro_t);
    return resu_t;
}
// Plot integrand as a function of p0 for several omegas, below 1, around 1 and above 1. 




// Plotting the integral 6.31 from Hatsuda 
double Pi_alpha(double p0, double omega,double p,double q,double T,double m,double eps){
    Int_gsl_adaptive int_t;
    double resu_t, erro_t;
    gsl_set_error_handler_off();
    funct h_lambda_t = [&] (double p){
        funct h_lambda1_t = [&] (double p0){
            int N_c = 3;
            double sign_p01 = (p0 + omega >= 0) ? 1.0 : -1.0;
            double sign_p02 = (p0 - omega >= 0) ? 1.0 : -1.0;
            complex<double> coeff1 = p * p * tanh(p0/2/T);
            complex<double> green1 = (p0 + omega)*(p0 + omega) - (p+q)*(p+q) - m*m + complex<double>(0,1) *eps ; 
            complex<double> green2 = p0*p0 - p*p - m*m + complex<double>(0,1) * eps;
            complex<double> green3 = (p0 - omega)*(p0 - omega) - (p - q)*(p - q) - m*m + complex<double>(0,1) *eps;

            return N_c 
        }
    }

}


