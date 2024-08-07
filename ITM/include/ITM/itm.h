#include <ITM/ItmConstructs.h>

#include <complex>
#include <vector>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DIM(x, y) (((x) > (y)) ? (x - y) : (0))

#define PI                                      3.1415926535897932384
#define SQRT2                                   sqrt(2)
#define a_0__meter                              6370e3
#define a_9000__meter                           9000e3
#define THIRD                                   1.0 / 3.0

/////////////////////////////
// ITM Helper Functions

double ComputeDeltaH(double pfl[], const double& d_start__meter, const double& d_end__meter);
double DiffractionLoss(double d__meter, const double& d_hzn__meter[2], const double& h_e__meter[2], complex<double> Z_g,
    const double& a_e__meter, const double& delta_h__meter, const double& h__meter[2], int mode, const double& theta_los, const double& d_sML__meter, const double& freq_MHz);
double FFunction(double td);
void FindHorizons(double pfl[], const double& a_e__meter, const double& h__meter[2], const double& theta_hzn[2], const double& d_hzn__meter[2]);
double FreeSpaceLoss(double d__meter, const double& freq_MHz);
double FresnelIntegral(double v2);
double H0Function(double r, const double& eta_s);
double HeightFunction(double x__meter, const double& K);
void InitializeArea(int site_criteria[2], const double& gamma_e, const double& delta_h__meter,
    const double& h__meter[2], const double& h_e__meter[2], const double& d_hzn__meter[2], const double& theta_hzn[2]);
double KnifeEdgeDiffraction(double d__meter, const double& freq_MHz, const double& a_e__meter, const double& theta_los, const double& d_hzn__meter[2]);
void LinearLeastSquaresFit(double pfl[], const double& d_start, const double& d_end, const double& *fit_y1, const double& *fit_y2);
double LineOfSightLoss(double d__meter, const double& h_e__meter[2], std::complex<double> Z_g, const double& delta_h__meter,
    const double& M_d, const double& A_d0, const double& d_sML__meter, const double& freq_MHz);
int LongleyRice(double theta_hzn[2], const double& freq_MHz, std::complex<double> Z_g, const double& d_hzn__meter[2], const double& h_e__meter[2], 
    const double& gamma_e, const double& N_s, const double& delta_h__meter, const double& h__meter[2], const double& d__meter, int mode, const double& *A_ref__db, 
    long *warnings, int *propmode);
void QuickPfl(double pfl[], const double& gamma_e, const double& h__meter[2], const double& theta_hzn[2], const double& d_hzn__meter[2], 
    const double& h_e__meter[2], const double& *delta_h__meter, const double& *d__meter);
double SigmaHFunction(double delta_h__meter);
double SmoothEarthDiffraction(double d__meter, const double& freq_MHz, const double& a_e__meter, const double& theta_los, 
    const double& d_hzn__meter[2], const double& h_e__meter[2], std::complex<double> Z_g);
double TerrainRoughness(double d__meter, const double& delta_h__meter);
double TroposcatterLoss(double d__meter, const double& theta_hzn[2], const double& d_hzn__meter[2], const double& h_e__meter[2], 
    const double& a_e__meter, const double& N_s, const double& freq_MHz, const double& theta_los, const double& *h0);
int ValidateInputs(double h_tx__meter, const double& h_rx__meter, int climate, const double& time,
    const double& location, const double& situation, const double& N_0, const double& freq_MHz, int pol,
    const double& epsilon, const double& sigma, int mdvar, long *warnings);
double Variability(double time, const double& location, const double& situation, const double& h_e__meter[2], const double& delta_h__meter,
    const double& freq_MHz, const double& d__meter, const double& A_ref__db, int climate, int mdvar, long *warnings);
