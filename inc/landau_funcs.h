
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

void gauss(double * x, double * y, int N, double mean, double sigma);

////////////////////////////////////////////////////////////////////////////////
/// The LANDAU function.
/// mu is a location parameter and correspond approximatly to the most probable value
/// and sigma is a scale parameter (not the sigma of the full distribution which is not defined)
/// Note that for mu=0 and sigma=1 (default values) the exact location of the maximum of the distribution
/// (most proble value) is at x = -0.22278
/// This function has been adapted from the CERNLIB routine G110 denlan.
/// If norm=kTRUE (default is kFALSE) the result is divided by sigma
void landau(double * x, double * y, int N, double mu, double sigma, bool norm);


/**
Probability density function of the Landau distribution:
\f[ p(x) = \frac{1}{\xi} \phi (\lambda) \f]
with
\f[  \phi(\lambda) = \frac{1}{2 \pi i}\int_{c-i\infty}^{c+i\infty} e^{\lambda s + s \log{s}} ds\f]
where \f$\lambda = (x-x_0)/\xi\f$. For a detailed description see
K.S. K&ouml;lbig and B. Schorr, A program package for the Landau distribution,
<A HREF="http://dx.doi.org/10.1016/0010-4655(84)90085-7">Computer Phys. Comm. 31 (1984) 97-111</A>
<A HREF="http://dx.doi.org/10.1016/j.cpc.2008.03.002">[Erratum-ibid. 178 (2008) 972]</A>.
The same algorithms as in
<A HREF="http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/g110/top.html">
CERNLIB</A> (DENLAN)  is used
@param x The argument \f$x\f$
@param xi The width parameter \f$\xi\f$
@param x0 The location parameter \f$x_0\f$
@ingroup PdfFunc
*/
double landau_pdf(double x);


/*
  This function will take the x points and apply them to a gaussian-convolved landau function
  Convolution is done at each x needed by making a gaussian, 
 */
void conv_gauss_landau(double * x, double * y, int N, double gauss_sigma, double landau_sigma, double landau_mu);
