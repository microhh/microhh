/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H

#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <algorithm>

// Prognostic = integrated in time with tendency, communicated between MPI tasks, and saved in restart files.
// Constant   = constant in time, communicated between MPI tasks, and saved in restart files.
// Diagnostic = recalculated every time step, not communicated or saved in restart files.
enum class Particle_kind {Prognostic, Constant, Diagnostic};


template<typename TF>
struct Particle_property
{
    Particle_kind kind;
    std::string long_name;
    std::string units;

    std::vector<TF> data;
    std::vector<TF> tend;   // Only used for prognostic properties.
};


template<typename TF>
class Particle_data
{
    public:
        void add(
                const std::string& name,
                const Particle_kind kind,
                const std::string& long_name,
                const std::string& units)
        {
            if (props.find(name) != props.end())
                throw std::runtime_error("Particle property \"" + name + "\" already exists");

            props.emplace(name, Particle_property<TF>{kind, long_name, units, {}, {}});
        }

        std::vector<TF>& get(const std::string& name)
        {
            return get_property(name).data;
        }

        std::vector<TF>& get_tend(const std::string& name)
        {
            auto& prop = get_property(name);

            if (prop.kind != Particle_kind::Prognostic)
                throw std::runtime_error("Particle property \"" + name + "\" is not prognostic and has no tendency");

            return prop.tend;
        }

        Particle_property<TF>& get_property(const std::string& name)
        {
            auto it = props.find(name);

            if (it == props.end())
                throw std::runtime_error("Particle property \"" + name + "\" does not exist");

            return it->second;
        }

        bool has_property(const std::string& name) const
        {
            return props.find(name) != props.end();
        }

        // Names of all properties of the requested kind(s), in alphabetical order.
        std::vector<std::string> get_names(std::initializer_list<Particle_kind> kinds) const
        {
            std::vector<std::string> names;

            for (const auto& [name, prop] : props)
                if (std::find(kinds.begin(), kinds.end(), prop.kind) != kinds.end())
                    names.push_back(name);

            return names;
        }

        int size() const
        {
            return uid.size();
        }

        // Resize all properties, and tendencies of prognostic properties.
        void resize(const int new_size)
        {
            uid.resize(new_size);

            for (auto& [name, prop] : props)
            {
                prop.data.resize(new_size);

                if (prop.kind == Particle_kind::Prognostic)
                    prop.tend.resize(new_size);
            }
        }

        std::vector<int> uid;

    private:
        // Use a `std::map`; the iteration order is identical on all MPI tasks, and
        // references to properties remain valid when new properties are added.
        std::map<std::string, Particle_property<TF>> props;
};
#endif
