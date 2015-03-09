/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
        master.start(argc, argv);

        // Print the current version of the model.
        master.print_message("Microhh git-hash: " GITHASH "\n");

        // Initialize the input class and read the input data from disk.
        Input input(&master);

        // Initialize the model class.
        Model model(&master, &input);

        // Initialize the master class.
        master.init(&input);

        // Initialize the model components.
        model.init();

        if (master.mode == "init")
        {
            // Initialize the allocated fields and save the data.
            model.save();
        }
        else if (master.mode == "run" || master.mode == "post")
        {
            // Initialize the allocated fields using data from disk.
            model.load();
        }

        // Print warnings for input variables that are unused.
        input.print_unused();

        // Free the memory taken by the input fields.
        input.clear();

        // Run the model.
        if (master.mode != "init")
            model.exec();
    }

    // Catch any exceptions and return 1.
    catch (...)
    {
        return 1;
    }

    // Return 0 in case of normal exit.
    return 0;
}
