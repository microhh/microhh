/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef MICROPHYS_SB06_BUDGET_H
#define MICROPHYS_SB06_BUDGET_H

#include <vector>
#include <string>
#include <map>

template<typename TF>
class Micro_budget {
    public:
        void create(const int kcells)
        {
            // All possible species.
            std::vector<std::string> species =
                    {"qv", "qc", "nc", "qi", "ni", "qr", "nr", "qs", "ns", "qg", "ng", "qh", "nh"};

            for (auto& specie : species)
                this->prev_tendencies.emplace(specie, std::vector<TF>(kcells));

            this->nlev = kcells;
        }

        void add_process(
                Stats<TF>& stats,
                const std::string& process,
                std::vector<std::string> species)
        {
            const std::string group_name = "sb06_budget";

            for (auto& specie : species)
            {
                const std::string full_name = specie + "_" + process;
                this->tendencies.emplace(full_name, std::vector<TF>(this->nlev));

                // Register in statistics.
                const std::string unit = full_name.at(0) == 'q' ? "kg kg-1 s" : "m-3 s";
                stats.add_prof(full_name, "Tendency " + full_name, unit, "z" , group_name);
            }
        }

        void set_stats(Stats<TF>& stats)
        {
            for (auto& prof : this->tendencies)
            {
                // Set time integrated tendency in stats.
                stats.set_prof(prof.first, prof.second);

                // Reset tendencies.
                std::fill(prof.second.begin(), prof.second.end(), TF(0));
            }
        }

        void set(
                const std::string& process,
                const std::string& specie,
                const TF tendency,
                const int k)
        {
            const std::string full_name = specie + "_" + process;
            this->tendencies.at(full_name)[k] += tendency - this->prev_tendencies.at(specie)[k];
            this->prev_tendencies.at(specie)[k] = tendency;
        }

        void reset_tendencies(const int k)
        {
            for (auto& prof : this->prev_tendencies)
                prof.second[k] = TF(0);
        }

    private:
        int nlev;
        std::map<std::string, std::vector<TF>> prev_tendencies;
        std::map<std::string, std::vector<TF>> tendencies;
};
#endif