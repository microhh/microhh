#ifndef OPTICAL_PROPS_H
#define OPTICAL_PROPS_H

#include <memory>
#include "Array.h"

template<typename TF>
class Optical_props
{
    public:
        Optical_props(
                const Array<TF,2>& band_lims_wvn,
                const Array<int,2>& band_lims_gpt);

        Optical_props(
                const Array<TF,2>& band_lims_wvn);

        virtual ~Optical_props() {};

        Optical_props(const Optical_props&) = default;

        Array<int,1> get_gpoint_bands() const { return this->gpt2band; }
        int get_nband() const { return this->band2gpt.dim(2); }
        int get_ngpt() const { return this->band2gpt.max(); }
        Array<int,2> get_band_lims_gpoint() const { return this->band2gpt; }
        Array<TF,2> get_band_lims_wavenumber() const { return this->band_lims_wvn; }

    private:
        Array<int,2> band2gpt;     // (begin g-point, end g-point) = band2gpt(2,band)
        Array<int,1> gpt2band;     // band = gpt2band(g-point)
        Array<TF,2> band_lims_wvn; // (upper and lower wavenumber by band) = band_lims_wvn(2,band)
};

// Base class for 1scl and 2str solvers fully implemented in header.
template<typename TF>
class Optical_props_arry : public Optical_props<TF>
{
    public:
        Optical_props_arry(const Optical_props<TF>& optical_props) :
            Optical_props<TF>(optical_props)
        {}

        virtual ~Optical_props_arry() {};

        virtual Array<TF,3>& get_tau() = 0;
        virtual Array<TF,3>& get_ssa() = 0;
        virtual Array<TF,3>& get_g  () = 0;

        virtual const Array<TF,3>& get_tau() const = 0;
        virtual const Array<TF,3>& get_ssa() const = 0;
        virtual const Array<TF,3>& get_g  () const = 0;

        virtual void set_subset(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
                const int col_s, const int col_e) = 0;

        virtual void get_subset(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
                const int col_s, const int col_e) = 0;

        virtual int get_ncol() const = 0;
        virtual int get_nlay() const = 0;
};

template<typename TF>
class Optical_props_1scl : public Optical_props_arry<TF>
{
    public:
        // Initializer constructor.
        Optical_props_1scl(
                const int ncol,
                const int nlay,
                const Optical_props<TF>& optical_props);

        void set_subset(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
                const int col_s, const int col_e);

        void get_subset(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
                const int col_s, const int col_e);

        int get_ncol() const { return tau.dim(1); }
        int get_nlay() const { return tau.dim(2); }

        Array<TF,3>& get_tau() { return tau; }
        Array<TF,3>& get_ssa() { throw std::runtime_error("ssa is not available in this class"); }
        Array<TF,3>& get_g  () { throw std::runtime_error("g is available in this class"); }

        const Array<TF,3>& get_tau() const { return tau; }
        const Array<TF,3>& get_ssa() const { throw std::runtime_error("ssa is not available in this class"); }
        const Array<TF,3>& get_g  () const { throw std::runtime_error("g is available in this class"); }

    private:
        Array<TF,3> tau;
};

template<typename TF>
class Optical_props_2str : public Optical_props_arry<TF>
{
    public:
        Optical_props_2str(
                const int ncol,
                const int nlay,
                const Optical_props<TF>& optical_props);

        void set_subset(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
                const int col_s, const int col_e);

        void get_subset(
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_sub,
                const int col_s, const int col_e);

        int get_ncol() const { return tau.dim(1); }
        int get_nlay() const { return tau.dim(2); }

        Array<TF,3>& get_tau() { return tau; }
        Array<TF,3>& get_ssa() { return ssa; }
        Array<TF,3>& get_g  () { return g; }

        const Array<TF,3>& get_tau() const { return tau; }
        const Array<TF,3>& get_ssa() const { return ssa; }
        const Array<TF,3>& get_g  () const { return g; }

    private:
        Array<TF,3> tau;
        Array<TF,3> ssa;
        Array<TF,3> g;
};

template<typename TF> void add_to(Optical_props_1scl<TF>& op_inout, const Optical_props_1scl<TF>& op_in);
template<typename TF> void add_to(Optical_props_2str<TF>& op_inout, const Optical_props_2str<TF>& op_in);
#endif
