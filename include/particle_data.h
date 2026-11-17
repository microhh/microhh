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
#include <cstring>

#ifdef USEMPI
#include <mpi.h>
#endif


template<typename T>
int adaptive_resize(
        std::vector<T>& v,
        const int new_size,
        const double reserve_ratio)
{
    const int size     = v.size();      // Used elements.
    const int capacity = v.capacity();  // Reserved size.

    // Keep track of increases and decreases in size.
    int status = 0;     // 0 = no change, 1 = increase, 2 = decrease.

    // Don't enlarge and shrink using the same ratio.
    const double shrink_ratio = 2. * reserve_ratio - 1.;

    if (new_size * shrink_ratio < capacity)
    {
        // Too large! Decrease capacity.
        status = 2;
        const int new_capacity = int(new_size * reserve_ratio);
        std::vector<T>(v.begin(), v.begin() + new_size).swap(v);
        v.reserve(new_capacity);

        //std::cout << "Too large! Size=" << size << " -> new size=" << new_size << ", new capacity=" << new_capacity << std::endl;
    }
    else if (new_size > capacity)
    {
        // Too small! Increase capacity.
        status = 1;
        const int new_capacity = int(new_size * reserve_ratio);
        v.reserve(new_capacity);

        //std::cout << "Too small! Size=" << size << " -> new size=" << new_size << ", new capacity=" << new_capacity << std::endl;
    }

    v.resize(new_size);

    return status;
}


// Particles are communicated between MPI tasks as records of bytes: [uid, field_0, field_1, ...],
// where the fields can be any list of particle properties and/or tendencies.
template<typename TF>
size_t get_record_size(const int n_fields)
{
    return sizeof(int) + n_fields * sizeof(TF);
}


template<typename TF>
void pack_record(
        char* buffer,
        const int n,
        const std::vector<int>& uid,
        const std::vector<std::vector<TF>*>& fields)
{
    std::memcpy(buffer, &uid[n], sizeof(int));
    buffer += sizeof(int);

    for (const std::vector<TF>* field : fields)
    {
        std::memcpy(buffer, &(*field)[n], sizeof(TF));
        buffer += sizeof(TF);
    }
}


template<typename TF>
void unpack_record(
        const char* buffer,
        const int n,
        std::vector<int>& uid,
        const std::vector<std::vector<TF>*>& fields)
{
    std::memcpy(&uid[n], buffer, sizeof(int));
    buffer += sizeof(int);

    for (std::vector<TF>* field : fields)
    {
        std::memcpy(&(*field)[n], buffer, sizeof(TF));
        buffer += sizeof(TF);
    }
}


// Pointers to a list of arrays, e.g. to pass temporary buffers as list of fields.
template<typename TF>
std::vector<std::vector<TF>*> get_pointers(std::vector<std::vector<TF>>& arrays)
{
    std::vector<std::vector<TF>*> pointers;

    for (std::vector<TF>& array : arrays)
        pointers.push_back(&array);

    return pointers;
}


#ifdef USEMPI
inline MPI_Datatype create_record_type(const size_t record_size)
{
    MPI_Datatype record_type;
    MPI_Type_contiguous(static_cast<int>(record_size), MPI_BYTE, &record_type);
    MPI_Type_commit(&record_type);
    return record_type;
}
#endif


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

        // Pointers to the data of properties `names`.
        std::vector<std::vector<TF>*> get_fields(const std::vector<std::string>& names)
        {
            std::vector<std::vector<TF>*> fields;

            for (const std::string& name : names)
                fields.push_back(&get(name));

            return fields;
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

        // Resize all arrays with `adaptive_resize()`.
        // Returns status of `uid`: 0 = no change, 1 = increase, 2 = decrease in capacity.
        int resize_adaptive(const int new_size, const TF reserve_ratio)
        {
            const int status = adaptive_resize(uid, new_size, reserve_ratio);

            for (auto& [name, prop] : props)
            {
                adaptive_resize(prop.data, new_size, reserve_ratio);

                if (prop.kind == Particle_kind::Prognostic)
                    adaptive_resize(prop.tend, new_size, reserve_ratio);
            }

            return status;
        }

        // Arrays that are communicated between MPI tasks: data of prognostic and
        // constant properties, and tendencies of prognostic properties.
        std::vector<std::vector<TF>*> get_comm_fields()
        {
            std::vector<std::vector<TF>*> fields;

            for (auto& [name, prop] : props)
            {
                if (prop.kind == Particle_kind::Diagnostic)
                    continue;

                fields.push_back(&prop.data);

                if (prop.kind == Particle_kind::Prognostic)
                    fields.push_back(&prop.tend);
            }

            return fields;
        }

        // Remove particles at the sorted indexes `leaving`, by moving particles from the end
        // of the arrays into the empty spots. The arrays are not resized.
        void compact(const std::vector<int>& leaving)
        {
            const std::vector<std::vector<TF>*> fields = get_comm_fields();

            // Create lookup table to check if an index is leaving.
            // Prevents using a find operation, which is more expensive.
            std::vector<char> is_leaving(uid.size(), 0);
            for (int idx : leaving)
                is_leaving[idx] = 1;

            int last_idx = uid.size()-1;

            for (int pos : leaving)
            {
                // If idx = [0,1,2,3,4,5] and `1` and `4` leave, the compacted vector is:
                //    idx = [0,5,2,3]
                // So if the leaving index is behind `last_idx`, there is no need to move it.
                if (pos > last_idx)
                    continue;

                // Find a valid index at the end of the vector.
                while (last_idx > pos && is_leaving[last_idx])
                    last_idx--;

                // Same as above; skip if already in empty block.
                if (pos >= last_idx)
                    continue;

                // Move data from end of vector to empty spot.
                uid[pos] = uid[last_idx];

                for (std::vector<TF>* field : fields)
                    (*field)[pos] = (*field)[last_idx];

                last_idx--;
            }
        }

        std::vector<int> uid;

    private:
        // Use a `std::map`; the iteration order is identical on all MPI tasks, and
        // references to properties remain valid when new properties are added.
        std::map<std::string, Particle_property<TF>> props;
};
#endif
