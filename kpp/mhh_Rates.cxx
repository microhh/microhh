/* Rate calculations */
template<typename TF>
TF EPR(const TF A1, const TF C1, const TF A2, const TF C2,
        const TF A3, const TF C3, const TF mmult, const TF ch2o,
        const TF TEMP)
{
    const TF K1 = A1 * exp(C1/TEMP);
    const TF K2 = A2 * exp(C2/TEMP) * mmult;
    const TF K3 = A3 * exp(C3/TEMP) * ch2o;
    const TF EPR_p = (K1 + K2) * (TF(1.0) + K3);

    return EPR_p;
}

template<typename TF>
TF ARR3(const TF A0, const TF B0, const TF TEMP)
{
    return A0 * exp(B0 / TEMP);
}

template<typename TF>
TF TROE_no2oh(
        const TF kzero, const TF mzero, const TF kinf,
        const TF fmulti, const TF MN2, const TF TEMP)
{
    const TF k0T  = (kzero * pow((TF(300)/TEMP), mzero)) * MN2;
    const TF kinfT = kinf;
    const TF znn = TF(0.75) - (TF(1.27) * log10(TF(0.41)));
    return (k0T * kinfT) / (k0T + kinfT) * pow(fmulti, (log10(TF(0.41)) / (TF(1) + pow((log10(k0T/kinfT))/znn, 2))));
}

template<typename TF>
TF TROE_cooh(
        const TF kzero, const TF mzero, double kinf,
        const TF minf, const TF k3, const TF c3, const TF k4,
        const TF c4, const TF fmulti, const TF MN2, const TF TEMP)
{
    const TF k0T  = (kzero * pow((TF(300)/TEMP), mzero)) * MN2;
    const TF kinfT = kinf * pow((TF(300)/TEMP), minf);
    const TF kval3 = k3 * pow((TF(300)/TEMP), c3);
    const TF kval4 = (k4 * pow((TF(300)/TEMP), c4)) / MN2;
    const TF kcooh = k0T / (TF(1) + k0T / kinfT) * pow(fmulti, TF(1) / (TF(1) + pow(log10(k0T/kinfT), 2)));
    return kcooh + (kval3 / (TF(1) + kval3 / kval4) * pow(fmulti, TF(1) / (TF(1) + pow(log10(kval3/kval4), 2))));
}

template<typename TF>
TF TROE_ifs(
        const TF kzero, const TF mzero, const TF kinf,
        const TF minf, const TF fmulti, const TF MN2, const TF TEMP)
{
    const TF k0T  = (kzero * pow((TEMP/TF(300)), mzero)) * MN2;
    const TF kinfT = kinf * pow((TEMP/TF(300)), minf);
    const TF znn_b = TF(0.75) - (TF(1.27) * log10(TF(0.35)));
    return (k0T * kinfT) / (k0T + kinfT) * pow(fmulti, log10(TF(0.35)) / (TF(1) + pow(log10(k0T/kinfT)/znn_b, 2)));
}

template<typename TF>
TF TROE_ifs2(
        const TF kzero, const TF mzero, const TF kinf,
        const TF minf, const TF fmulti, const TF MN2, const TF c1,
        const TF c2, const TF TEMP)
{
    const TF k0T  = (kzero * pow((TEMP/TF(300)), mzero)) * MN2 * exp(c1 / TEMP);
    const TF kinfT = kinf * pow((TEMP/TF(300)), minf) * exp(c2 / TEMP);
    const TF znn_b = TF(0.75) - (TF(1.27) * log10(TF(0.35)));
    return (k0T * kinfT) / (k0T + kinfT) * pow(fmulti, log10(TF(0.35)) / (TF(1) + pow(log10(k0T/kinfT)/znn_b, 2)));
}

template<typename TF>
TF k3rd_iupac(
        const TF kzero, const TF mzero, const TF kinf,
        const TF minf, const TF fmulti, const TF MN2,
        const TF NM, const TF TEMP)
{
    const TF k0T  = (kzero * pow((TEMP/TF(300)), mzero)) * MN2;
    const TF kinfT = kinf * pow((TEMP/TF(300)), minf);
    return (k0T / (TF(1) + (k0T/kinfT))) * pow( fmulti, TF(1)/ (TF(1) + pow((log10(k0T/kinfT)/NM), 2)));
}

template<typename TF>
TF usr_O3_hv_H2O(const TF TEMP, const TF C_M, const TF C_H2O, const TF J_O1D)
{
    const TF KH2O = (TF(1.63e-10) * C_H2O * exp(TF(60)/TEMP)  )  ;
    const TF KN2  = (TF(2.15e-11) * exp(TF(110.0)/TEMP) * TF(0.79)*C_M) ;
    const TF KO2  = (TF(3.30e-11) * exp(TF(55.0) /TEMP) * TF(0.21)*C_M) ;
    return (KH2O *J_O1D) / (KH2O + KN2 + KO2);
}

template<typename TF>
TF RK28(
        const TF k0a, const TF k0ea, const TF k2a,
        const TF k2ea, const TF k3a, const TF k3ea,
        const TF MN2, const TF TEMP)
{
    const TF k0 = k0a * exp(k0ea/TEMP);
    const TF k2 = (k2a * exp(k2ea/TEMP)) * MN2;
    const TF k3 = k3a * exp(k3ea/TEMP);
    return k0 + k3 * k2 / (k3 + k2);
}

template<typename TF>
TF EPR_ev(
        const TF A1, const TF C1,
        const TF A2, const TF C2,
        const TF mmult, const TF TEMP)
{
    const TF K1 = A1 * exp(C1/TEMP);
    const TF K2 = A2 * exp(C2/TEMP) * mmult;
    const TF EPR_p = (K1 + K2);
    return EPR_p;
}

template<typename TF>
TF kohch4(const TF A, const TF b, const TF C, const TF TEMP)
{
    return A * pow(TEMP, b) * exp(C/TEMP);
}

template<typename TF>
void isop_stat(
        const TF* tisop, const TF* isop,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int kstart, const int kend,
        const int jj, const int kk)
{
    for (int k=kstart; k<kend; ++k)
    {
        TF iso_mean = TF(0.0);
        TF iso_tend = TF(0.0);
        int cnt = 0;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                iso_mean += isop[ijk];
                iso_tend += tisop[ijk];
                cnt += 1;
            }
        //printf("%i  %12.3e ppb  %12.3e ppb/hour \n",k,iso_mean/cnt,iso_tend*3600.0/cnt);
    }
}
