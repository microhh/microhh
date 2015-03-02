namespace most
{
  //
  // GRADIENT FUNCTIONS
  //
  inline double phim_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return std::pow(1. + 3.6*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }

  inline double phim_stable(const double zeta)
  {
    // Hogstrom, 1988
    return 1. + 4.8*zeta;
  }

  inline double phim(const double zeta)
  {
    return (zeta <= 0.) ? phim_unstable(zeta) : phim_stable(zeta);
  }

  inline double phih_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }

  inline double phih_stable(const double zeta)
  {
    // Hogstrom, 1988
    return 1. + 7.8*zeta;
  }

  inline double phih(const double zeta)
  {
    return (zeta <= 0.) ? phih_unstable(zeta) : phih_stable(zeta);
  }

  //
  // INTEGRATED FUNCTIONS
  //
  inline double psim_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return 3.*std::log( ( 1. + 1./phim_unstable(zeta) ) / 2.);
  }

  inline double psim_stable(const double zeta)
  {
    // Hogstrom, 1988
    return -4.8*zeta;
  }

  inline double psih_unstable(const double zeta)
  {
    // Wilson, 2001 functions, see Wyngaard, page 222.
    return 3. * std::log( ( 1. + 1. / phih_unstable(zeta) ) / 2.);
  }

  inline double psih_stable(const double zeta)
  {
    // Hogstrom, 1988
    return -7.8*zeta;
  }

  inline double fm(const double zsl, const double z0m, const double L)
  {
    return (L <= 0.)
      ? constants::kappa / (std::log(zsl/z0m) - psim_unstable(zsl/L) + psim_unstable(z0m/L))
      : constants::kappa / (std::log(zsl/z0m) - psim_stable  (zsl/L) + psim_stable  (z0m/L));
  }

  inline double fh(const double zsl, const double z0h, const double L)
  {
    return (L <= 0.)
      ? constants::kappa / (std::log(zsl/z0h) - psih_unstable(zsl/L) + psih_unstable(z0h/L))
      : constants::kappa / (std::log(zsl/z0h) - psih_stable  (zsl/L) + psih_stable  (z0h/L));
  }
}
