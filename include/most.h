namespace most
{
  inline double psim(const double zeta)
  {
    double psim;
    if(zeta <= 0.)
    {
      // Wilson, 2001 functions, see Wyngaard, page 222.
      const double x = std::pow(1. + std::pow(3.6 * std::abs(zeta),2./3.), -0.5);
      psim = 3.*std::log( (1. + 1./x) / 2.);
    }
    else
    {
      // Hogstrom, 1988
      psim = -4.8*zeta;
    }
    return psim;
  }
  
  inline double psih(const double zeta)
  {
    double psih;
    if(zeta <= 0.)
    {
      // Wilson, 2001 functions, see Wyngaard, page 222.
      const double x = std::pow(1. + std::pow(7.9*std::abs(zeta), (2./3.)), -0.5);
      psih = 3. * std::log( (1. + 1. / x) / 2.);
    }
    else
    {
      // Hogstrom, 1988
      psih  = -7.8*zeta;
    }
    return psih;
  }
  
  inline double phim(const double zeta)
  {
    double phim;
    if(zeta <= 0.)
    {
      // Wilson, 2001 functions, see Wyngaard, page 222.
      phim = std::pow(1. + 3.6*std::pow(std::abs(zeta), 2./3.), -1./2.);
    }
    else
      // Hogstrom, 1988
      phim = 1. + 4.8*zeta;
  
    return phim;
  }
  
  inline double phih(const double zeta)
  {
    double phih;
    if(zeta <= 0.)
    {
      // Wilson, 2001 functions, see Wyngaard, page 222.
      phih = std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
    }
    else
      // Hogstrom, 1988
      phih = 1. + 7.8*zeta;
  
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
