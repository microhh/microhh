#include <cmath>
#include <iostream>
#include <fstream>

// Compile this file to produce reference data for the test case.

// Input parameters
namespace
{
    const double L = 5.12;
    const int nz = 1025;       // Number of grid points in z direction. So 1024 intervals.
    const int nx = 513;        // Number of grid points in x direction. So 512 intervals.
    const double nu = 0.001;   // Kinematic viscosity.
    const double alp = 0.001;  // Thermal diffusivity, alpha.
    const double BV = 0.02;    // Brunt-Vaisala frequency N.
    const double bmax = 1e-5;  // Max surface buoyancy.
    const int nmax = 50000;    // Number of terms in the Fourier series. Paper user 50000.

    // Derived variables
    const double dx = L / (nx-1); // Horizontal grid spacing in meters.
    const double dz = dx;         // Vertical grid spacing identical to horizontal.

    const double pi = 4.*std::atan(1.);

    const int ntot = nx*nz;
    const double dx2 = dx*dx;
    const double dz2 = dz*dz;

    const double BV43 = std::pow(BV, 4./3.);
    const double bmin = -bmax; // Min surface buoyancy; set to bmin = -bmax.

    const double arg = 2.*pi/3.;
    const double cs = std::cos(arg);
    const double sn = std::sin(arg);
    const double denom = std::pow(nu*alp, 1./3.);
}

int main()
{
    // Allocate and initialize arrays.
    double* z = new double[nz];
    double* x = new double[nx];

    double* b   = new double[ntot];
    double* psi = new double[ntot];
    double* w   = new double[ntot];
    double* eta = new double[ntot];
    double* u   = new double[ntot];

    for (int k=0; k<nz; ++k)
        z[k] = k*dz;

    for (int i=0; i<nx; ++i)
        x[i] = i*dx;

    // Set all arrays to zero.
    for (int n=0; n<ntot; ++n)
    {
        b   [n] = 0.;
        psi [n] = 0.;
        w   [n] = 0.;
        eta [n] = 0.;
        u   [n] = 0.;
    }

    // Calculate the solution on the grid.
    for (int n=1; n<nmax+1; ++n) // Start Fourier loop.
    {
        if (n % 100 == 0)
            std::cout << "n = " << n << std::endl;
        const double kappa = n*pi/L; // Define the wave number kappa.
        const double kappa2 = kappa*kappa;

        const double bn = 2. * (bmax*(1.-std::cos(n*pi/2.)) + bmin*(std::cos(n*pi/2.)-std::cos(n*pi)))
                             / (n*pi);

        // Compute M0, r, and phi   
        const double numer = std::pow(BV*kappa, 2./3.);
        const double M0 = - std::sqrt(kappa2 + numer/denom);
        const double r = std::sqrt( std::pow(kappa2 + numer*cs/denom, 2) + std::pow(numer*sn/denom, 2) );
        const double phi = std::acos( (kappa2 + numer*cs/denom)/r );

        const double mu = M0 / std::sqrt(r);

        const double denom2 = mu + 2.*std::cos(pi/3. + phi/2.);

        for (int k=0; k<nz; ++k)
        {
            const double Zs = z[k]*std::sqrt(r)*std::sin(phi/2.);
            const double Zc = z[k]*std::sqrt(r)*std::cos(phi/2.);

            const double arg1 = - Zc;
            const double arg2 = M0*z[k];

            const double exp1 = arg1 < -200. ? 0. : std::exp(arg1);
            const double exp2 = arg2 < -200. ? 0. : std::exp(arg2);

            const double numer_b = exp1 * ( mu*std::cos(Zs+pi/6.)
                                          + std::cos(Zs+pi/6.+phi/2.)) 
                                 - exp2 * std::sin(phi/2.);
            const double fac_b = 2.*bn/sqrt(3.);

            const double numer_psi = exp1 * ( mu*std::sin(Zs) 
                                            + std::sin(Zs+phi/2.))
                                   - exp2 * std::sin(phi/2.);

            const double fac_psi = (2.*bn*std::pow(alp, 2./3.))
                                 / (std::sqrt(3.)*BV43*std::pow(kappa*nu, 1./3.));

            const double fac_w = (2.*bn*std::pow(alp*kappa, 2./3.))
                               / (sqrt(3.)*BV43*std::pow(nu, 1./3.));
   
            const double numer_eta = exp1 * (mu*std::cos(Zs-pi/6.) + std::cos(Zs-pi/6.+phi/2.))
                                   + exp2 * std::sin(phi/2.);
            const double fac_eta = - (2.*bn*std::pow(alp*kappa, 1./3.))
                                 / (std::sqrt(3.)*std::pow(BV*nu, 2./3.));

            for (int i=0; i<nx; ++i)
            {
                const int ik = i + k*nx;
                const double fac_sin = std::sin(kappa*x[i])/denom2;
                const double fac_cos = std::cos(kappa*x[i])/denom2;
                b  [ik] += fac_b   * numer_b   * fac_sin;
                psi[ik] += fac_psi * numer_psi * fac_cos;
                w  [ik] += fac_w   * numer_psi * fac_sin;
                eta[ik] += fac_eta * numer_eta * fac_cos;
            } // i
        } // k
    } // n

    // Finite difference psi to get u. Skip top and bottom points.
    for (int k=1; k<nz-1; ++k)
    {
        for (int i=0; i<nx; ++i)
        {
            const int ik = i + k*nx;
            u[ik] = (psi[ik+nx] - psi[ik-nx]) / (2.*dz);
        }
    }

    // Verify that surface buoyancy computed from the Fourier series is really a square wave.
    // for (int i=0; i<nx; ++i)
    //     std::cout << "i, x, bs = " << i << ", " << x[i] << ", " << b[i] << std::endl;

    // Write the output data.
    std::ofstream uFile("u.out", std::ios::out | std::ios::binary);
    char* uchar = reinterpret_cast<char*>(u);
    uFile.write(uchar, ntot*sizeof(double));
    uFile.close();

    std::ofstream wFile("w.out", std::ios::out | std::ios::binary);
    char* wchar = reinterpret_cast<char*>(w);
    wFile.write(wchar, ntot*sizeof(double));
    wFile.close();

    std::ofstream bFile("b.out", std::ios::out | std::ios::binary);
    char* bchar = reinterpret_cast<char*>(b);
    bFile.write(bchar, ntot*sizeof(double));
    bFile.close();

    std::ofstream psiFile("psi.out", std::ios::out | std::ios::binary);
    char* psichar = reinterpret_cast<char*>(psi);
    psiFile.write(psichar, ntot*sizeof(double));
    psiFile.close();

    std::ofstream etaFile("eta.out", std::ios::out | std::ios::binary);
    char* etachar = reinterpret_cast<char*>(eta);
    etaFile.write(etachar, ntot*sizeof(double));
    etaFile.close();

    // Delete arrays.
    delete[] z;
    delete[] x;

    delete[] b;
    delete[] psi;
    delete[] w;
    delete[] eta;
    delete[] u;

    return 0;
};
