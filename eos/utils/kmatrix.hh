#ifndef EOS_GUARD_EOS_UTILS_KMATRIX_HH
#define EOS_GUARD_EOS_UTILS_KMATRIX_HH 1


#include <eos/utils/complex.hh>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <array>

namespace eos

{
  template <unsigned nchannels_, unsigned nresonances_>
  class KMatrix
  {

  public:

    struct Channel;
    struct Resonance;

    std::array<KMatrix::Channel, nchannels_> _channels;
    std::array<KMatrix::Resonance, nresonances_> _resonances;
    const std::string & _prefix;

    // Khat contains the normalized K-matrix entries
    // cf. S. U. Chung et al., Annalen Phys. 4, 404 (1995).
    gsl_matrix_complex * _Khat;
    // That is defined as T = (1-i*rho*K)^(-1)*K
    // (note that the two factors commute in the product)
    gsl_matrix_complex * _That;

    // Two temporary matrices and a permutation are needed to evaluate That
    gsl_matrix_complex * _tmp_1;
    gsl_matrix_complex * _tmp_2;
    gsl_permutation * _perm;

    // Constructor
    // KMatrix(std::initializer_list<KMatrix::Channel> channels,
    //         std::initializer_list<KMatrix::Resonance> resonances,
    //         const std::string & prefix);
    KMatrix(std::array<KMatrix::Channel, nchannels_> channels,
            std::array<KMatrix::Resonance, nresonances_> resonances,
            const std::string & prefix);

    
    // Destuctor
    ~KMatrix();

    // Return rowindex^th row of the T matrix defined as T = (1-i*rho*K)^(-1)*K
    std::array<complex<double>, nchannels_> tmatrix_row(unsigned rowindex,
                                                        const double s) const;

  };

  template <unsigned nchannels_, unsigned nresonances_>
  struct KMatrix<nchannels_, nresonances_>::Channel
  {
    std::string name;

    //Masses of the two final state particles
    double m1;
    double m2;

    // Normalized couplings to the resonances
    virtual std::array<double, nresonances_> g0() const = 0;

    // Phase space factor
    virtual double beta(const double & s) const = 0;
    // Analytic continuation of the phase space factor
    virtual complex<double> rho(const double & s) const = 0;
  };

  template <unsigned nchannels_, unsigned nresonances_>
  struct KMatrix<nchannels_, nresonances_>::Resonance
  {
    std::string name;

    //Mass of the resonance
    double m;
  };

}

#endif
