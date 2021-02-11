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

#include <iostream>
#include "master.h"
#include "model.h"

int main(int argc, char *argv[])
{
    // Initialize the master class, it cannot fail.
    Master master;
    try
    {
        // Start up the master class and the Message Passing Interface.
        master.start();

        // Print the current version of the model.
        master.print_message("Microhh git-hash: " GITHASH "\n");

        // Initialize the model in precision.
        #ifdef FLOAT_SINGLE
        master.print_message("Precision: Single (32-bits floats)\n");
        Model<float> model(master, argc, argv);
        #else
        master.print_message("Precision: Double (64-bits floats)\n");
        Model<double> model(master, argc, argv);
        #endif

        // Initialize the model components.
        model.init();

        // Load or save the data depending on the run mode.
        model.load_or_save();

        // Run the model.
        model.exec();
    }

    // Catch any exceptions and return 1.
    catch (const std::exception& e)
    {
        master.print_message("EXCEPTION: %s\n", e.what());
        return 1;
    }
    catch (...)
    {
        master.print_message("UNHANDLED EXCEPTION!\n");
        return 1;
    }

    // Return 0 in case of normal exit.
    return 0;
}
