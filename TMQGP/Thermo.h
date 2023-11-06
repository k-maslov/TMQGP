#ifndef THERMO_H_
#define THERMO_H_
#include "Interpolator.h"
#include "IntGSL.h"
#include <complex>

double OmQ_F_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG);
double OmQ_F(double T, Interpolator2D & iImG, Interpolator2D & iReG);
double OmQ_B_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG);
double OmQ_B(double T, Interpolator2D & iImG, Interpolator2D & iReG);
double delta(double om, double q, Interpolator2D & iImG, Interpolator2D & iReG);

double OmS_F_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);

double OmS_F(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);

double OmS_B_om_int(double q, double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);

double OmS_B(double T, Interpolator2D & iImG, Interpolator2D & iReG,
                                        Interpolator2D & iImS, Interpolator2D & iReS);

#endif