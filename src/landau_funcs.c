#include "landau_funcs.h"

void gauss(double * x, double * y, int N, double mean, double sigma) {
  static const float inv_sqrt_2pi = 0.3989422804014327;
  int n;
  for (n = 0; n < N; n++) {
    float a = (x[n] - mean) / sigma;
    y[n] = inv_sqrt_2pi / sigma * exp(-0.5f * a * a);
  }
  return;
}

////////////////////////////////////////////////////////////////////////////////
/// The LANDAU function.
/// mu is a location parameter and correspond approximatly to the most probable value
/// and sigma is a scale parameter (not the sigma of the full distribution which is not defined)
/// Note that for mu=0 and sigma=1 (default values) the exact location of the maximum of the distribution
/// (most proble value) is at x = -0.22278
/// This function has been adapted from the CERNLIB routine G110 denlan.
/// If norm=kTRUE (default is kFALSE) the result is divided by sigma

void landau(double * x, double * y, int N, double mu, double sigma, bool norm)
{
  int n;
  for (n = 0; n < N; n++) {
    if (sigma <= 0) {
      y[n] = 0;
    }
    double den = landau_pdf( (x[n] - mu) / sigma );
    if (!norm) {
      y[n] = den;
    }
    else {
      y[n] = den / sigma;
    }
  }
  return;
}

double landau_pdf(double x) {
  double xi = 1.0;
  double x0 = 0.0;

  // LANDAU pdf : algorithm from CERNLIB G110 denlan
  // same algorithm is used in GSL

  static double p1[5] = {0.4259894875, -0.1249762550, 0.03984243700, -0.006298287635,   0.001511162253};
  static double q1[5] = {1.0         , -0.3388260629, 0.09594393323, -0.01608042283,    0.003778942063};

  static double p2[5] = {0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411,   0.0001283617211};
  static double q2[5] = {1.0         , 0.7428795082, 0.3153932961,   0.06694219548,    0.008790609714};

  static double p3[5] = {0.1788544503, 0.09359161662, 0.006325387654, 0.00006611667319, -0.000002031049101};
  static double q3[5] = {1.0         , 0.6097809921, 0.2560616665,   0.04746722384,    0.006957301675};

  static double p4[5] = {0.9874054407, 118.6723273,  849.2794360,   -743.7792444,      427.0262186};
  static double q4[5] = {1.0         , 106.8615961,  337.6496214,    2016.712389,      1597.063511};

  static double p5[5] = {1.003675074,  167.5702434,  4789.711289,    21217.86767,     -22324.94910};
  static double q5[5] = {1.0         , 156.9424537,  3745.310488,    9834.698876,      66924.28357};

  static double p6[5] = {1.000827619,  664.9143136,  62972.92665,    475554.6998,     -5743609.109};
  static double q6[5] = {1.0         , 651.4101098,  56974.73333,    165917.4725,     -2815759.939};

  static double a1[3] = {0.04166666667, -0.01996527778, 0.02709538966};

  static double a2[2] = { -1.845568670, -4.284640743};

  if (xi <= 0) return 0;
  double v = (x - x0) / xi;
  double u, ue, us, denlan;
  if (v < -5.5) {
    u   = exp(v + 1.0);
    if (u < 1e-10) return 0.0;
    ue  = exp(-1 / u);
    us  = sqrt(u);
    denlan = 0.3989422803 * (ue / us) * (1 + (a1[0] + (a1[1] + a1[2] * u) * u) * u);
  } else if (v < -1) {
    u   = exp(-v - 1);
    denlan = exp(-u) * sqrt(u) *
             (p1[0] + (p1[1] + (p1[2] + (p1[3] + p1[4] * v) * v) * v) * v) /
             (q1[0] + (q1[1] + (q1[2] + (q1[3] + q1[4] * v) * v) * v) * v);
  } else if (v < 1) {
    denlan = (p2[0] + (p2[1] + (p2[2] + (p2[3] + p2[4] * v) * v) * v) * v) /
             (q2[0] + (q2[1] + (q2[2] + (q2[3] + q2[4] * v) * v) * v) * v);
  } else if (v < 5) {
    denlan = (p3[0] + (p3[1] + (p3[2] + (p3[3] + p3[4] * v) * v) * v) * v) /
             (q3[0] + (q3[1] + (q3[2] + (q3[3] + q3[4] * v) * v) * v) * v);
  } else if (v < 12) {
    u   = 1 / v;
    denlan = u * u * (p4[0] + (p4[1] + (p4[2] + (p4[3] + p4[4] * u) * u) * u) * u) /
             (q4[0] + (q4[1] + (q4[2] + (q4[3] + q4[4] * u) * u) * u) * u);
  } else if (v < 50) {
    u   = 1 / v;
    denlan = u * u * (p5[0] + (p5[1] + (p5[2] + (p5[3] + p5[4] * u) * u) * u) * u) /
             (q5[0] + (q5[1] + (q5[2] + (q5[3] + q5[4] * u) * u) * u) * u);
  } else if (v < 300) {
    u   = 1 / v;
    denlan = u * u * (p6[0] + (p6[1] + (p6[2] + (p6[3] + p6[4] * u) * u) * u) * u) /
             (q6[0] + (q6[1] + (q6[2] + (q6[3] + q6[4] * u) * u) * u) * u);
  } else {
    u   = 1 / (v - v * log(v) / (v + 1));
    denlan = u * u * (1 + (a2[0] + a2[1] * u) * u);
  }
  return denlan / xi;

}

void conv_gauss_landau(double * x, double * y, int N, double gauss_sigma, double landau_sigma, double landau_mu) {

  // Need to make a gaussian function centered at each x value, then convolve it with a landau function

  // Make a set of x values that covers the same range of x and at least the granularity, with a minimum
  // of 100 points


  int n_points = 200;
  double min;

  // Allocate the memory here at the beginning, and free at the very end.
  double * _internal_x = (double *) malloc(n_points * sizeof(double));

  double * _internal_gaussian = (double *) malloc(n_points * sizeof(double));
  double * _internal_landau   = (double *) malloc(n_points * sizeof(double));

  // This block of code gets the convolution of a gaussian and a landau at a single point x

  double step_size = (10.0 * gauss_sigma) / n_points;

  int n;
  for (n = 0; n < N; n++) {
    // Need to make a gaussian and landau distribution here.
    // The gaussian is centered on x[n], with width gauss_sigma, over _internal_x
    // The landua is over the same range

    min = x[n] - 5 * gauss_sigma;
    // Fill the points along that axis
    int i;
    for (i = 0; i < n_points; i ++) {
      _internal_x[i] = min + step_size * i;
    }

    // Fill the arrays using the other functions:
    gauss(_internal_x, _internal_gaussian, n_points, x[n], gauss_sigma);
    landau(_internal_x, _internal_landau, n_points, landau_mu, landau_sigma, true);

    // Now convolve the two.
    // In this case, since the gaussian is not summed to 1, multiply by the step size at each point
    y[n] = 0.0;
    for (i = 0; i < n_points; i ++){
      y[n] += _internal_gaussian[i] * _internal_landau[i] * step_size;
    }    

  }

  free(_internal_x);
  free(_internal_landau);
  free(_internal_gaussian);
  return;

}


