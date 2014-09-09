/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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
  // initialize the master class, it cannot fail
  cmaster master;
  try
  {
    // start up the master class
    master.startup(argc, argv);

    // print the current status of the model
    master.printMessage("Microhh git-hash: " GITHASH "\n");

    // create the instances of the objects
    cinput  input (&master);
    cmodel  model (&master, &input);

    // read the input data
    if(input.readinput())
      return 1;

    if(master.readinifile(&input))
      return 1;
    if(model.readinifile())
      return 1;

    // init the mpi 
    if(master.init())
      return 1;
    if(model.init())
      return 1;

    if(master.mode == "init")
    {
      if(model.create())
        return 1;

      // save the data
      if(model.save())
        return 1;
    }
    else
    {
      // load the data
      if(model.load())
        return 1;
    }

    // check unused input
    input.printUnused(); 
    // free the memory of the input
    input.clear();

    // run the model
    if(master.mode != "init")
      if(model.exec())
        return 1;
  }

  // catch any exceptions and return 1
  catch (...)
  {
    return 1;
  }

  // return 0 in case the model exits properly
  return 0;
}

