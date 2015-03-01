namespace most
{
  inline double psim(const double zeta)
  {
    double psim;
    if(zeta <= 0.)
    {
      // Businger-Dyer functions
      // const double x     = (1. - 16. * zeta) ** (0.25)
      // psim  = 3.14159265 / 2. - 2. * arctan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
      // Wilson functions
      const double x = std::pow(1. + std::pow(3.6 * std::abs(zeta),2./3.), -0.5);
      psim = 3.*std::log( (1. + 1./x) / 2.);
    }
    else
    {
      psim = -2./3.*(zeta - 5./0.35) * std::exp(-0.35 * zeta) - zeta - (10./3.) / 0.35;
    }
    return psim;
  }
  
  inline double psih(const double zeta)
  {
    double psih;
    if(zeta <= 0.)
    {
      // Businger-Dyer functions
      // const double x     = (1. - 16. * zeta) ** (0.25)
      // psih  = 2. * log( (1. + x ** 2.) / 2. )
      // Wilson functions
      const double x = std::pow(1. + std::pow(7.9*std::abs(zeta), (2./3.)), -0.5);
      psih = 3. * std::log( (1. + 1. / x) / 2.);
    }
    else
    {
      psih  = (-2./3.) * (zeta-5./0.35) * std::exp(-0.35*zeta) - std::pow(1. + (2./3.) * zeta, 1.5) - (10./3.) / 0.35 + 1.;
    }
    return psih;
  }
  
  inline double phim(const double zeta)
  {
    double phim;
    if(zeta <= 0.)
    {
      // Businger-Dyer functions
      // phim  = (1. - 16. * zeta) ** (-0.25)
      // Wilson functions
      phim = std::pow(1. + 3.6*std::pow(std::abs(zeta), 2./3.), -1./2.);
    }
    else
      phim = 1. + 5.*zeta;
  
    return phim;
  }
  
  inline double phih(const double zeta)
  {
    double phih;
    if(zeta <= 0.)
    {
      // Businger-Dyer functions
      // phih  = (1. - 16. * zeta) ** (-0.5)
      // Wilson functions
      phih = std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
    }
    else
      phih = 1. + 5.*zeta;
  
    return phih;
  }

  inline double fm(const double zsl, const double z0m, const double L)
  {
    return constants::kappa / (std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L));
  }
  
  inline double fh(const double zsl, const double z0h, const double L)
  {
    return constants::kappa / (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L));
  }
} 
