#ifndef FLUXES_H
#define FLUXES_H

#include <memory>

// Forward declarations.
template<typename TF, int> class Array;
template<typename TF> class Optical_props_arry;

template<typename TF>
class Fluxes
{
    public:
        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1) = 0;

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const Array<TF,3>& gpt_flux_dn_dir,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1) = 0;
};

template<typename TF>
class Fluxes_broadband : public Fluxes<TF>
{
    public:
        Fluxes_broadband(const int ncol, const int nlev);
        virtual ~Fluxes_broadband() {};

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1);

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const Array<TF,3>& gpt_flux_dn_dir,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1);

        Array<TF,2>& get_flux_up    () { return flux_up;     }
        Array<TF,2>& get_flux_dn    () { return flux_dn;     }
        Array<TF,2>& get_flux_dn_dir() { return flux_dn_dir; }
        Array<TF,2>& get_flux_net   () { return flux_net;    }

        virtual Array<TF,3>& get_bnd_flux_up    () { throw std::runtime_error("Band fluxes are not available"); }
        virtual Array<TF,3>& get_bnd_flux_dn    () { throw std::runtime_error("Band fluxes are not available"); }
        virtual Array<TF,3>& get_bnd_flux_dn_dir() { throw std::runtime_error("Band fluxes are not available"); }
        virtual Array<TF,3>& get_bnd_flux_net   () { throw std::runtime_error("Band fluxes are not available"); }

    private:
        Array<TF,2> flux_up;
        Array<TF,2> flux_dn;
        Array<TF,2> flux_dn_dir;
        Array<TF,2> flux_net;
};

template<typename TF>
class Fluxes_byband : public Fluxes_broadband<TF>
{
    public:
        Fluxes_byband(const int ncol, const int nlev, const int nbnd);
        virtual ~Fluxes_byband() {};

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1);

        virtual void reduce(
                const Array<TF,3>& gpt_flux_up,
                const Array<TF,3>& gpt_flux_dn,
                const Array<TF,3>& gpt_flux_dn_dir,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props,
                const int top_at_1);

        Array<TF,3>& get_bnd_flux_up    () { return bnd_flux_up;     }
        Array<TF,3>& get_bnd_flux_dn    () { return bnd_flux_dn;     }
        Array<TF,3>& get_bnd_flux_dn_dir() { return bnd_flux_dn_dir; }
        Array<TF,3>& get_bnd_flux_net   () { return bnd_flux_net;    }

    private:
        Array<TF,3> bnd_flux_up;
        Array<TF,3> bnd_flux_dn;
        Array<TF,3> bnd_flux_dn_dir;
        Array<TF,3> bnd_flux_net;
};
#endif
