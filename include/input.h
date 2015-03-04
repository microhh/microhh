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

#ifndef INPUT
#define INPUT

#include <map>
#include <string>
#include <vector>

// Forward declaration to avoid circular dependency.
class Master;

typedef std::map<std::string, std::vector<double> > DataMap;

class Input
{
  public:
    Input(Master *);
    ~Input();

    void clear();

    // Item retrieval functions
    int getItem(int *        , std::string, std::string, std::string);
    int getItem(int *        , std::string, std::string, std::string, int);
    int getItem(double *     , std::string, std::string, std::string);
    int getItem(double *     , std::string, std::string, std::string, double);
    int getItem(bool *       , std::string, std::string, std::string);
    int getItem(bool *       , std::string, std::string, std::string, bool);
    int getItem(std::string *, std::string, std::string, std::string);
    int getItem(std::string *, std::string, std::string, std::string, std::string);

    // List retrieval functions
    int getList(std::vector<int> *        , std::string, std::string, std::string);
    int getList(std::vector<double> *     , std::string, std::string, std::string);
    int getList(std::vector<std::string> *, std::string, std::string, std::string);

    int get_prof(double *, std::string, int size);
    int getTime(double **, std::vector<double> *, std::string);
    int getTimeProf(double **, std::vector<double> *, std::string, int);

    void printUnused();
    void flagUsed(std::string, std::string);

  private:
    Master *master;

    int readIniFile();
    int readDataFile(DataMap *, std::string, bool);

    template <class valuetype>
    int parseItem(valuetype *, std::string, std::string, std::string, bool, valuetype);

    template <class valuetype>
    int parseList(std::vector<valuetype> *, std::string, std::string, std::string);

    int checkItemExists(std::string, std::string, std::string el="default");
    int checkItem(int *        , std::string, std::string, std::string el="default");
    int checkItem(double *     , std::string, std::string, std::string el="default");
    int checkItem(bool *       , std::string, std::string, std::string el="default");
    int checkItem(std::string *, std::string, std::string, std::string el="default");

    // list retrieval
    int checkList(std::vector<int> *        , std::string, std::string, std::string el="default");
    int checkList(std::vector<double> *     , std::string, std::string, std::string el="default");
    int checkList(std::vector<std::string> *, std::string, std::string, std::string el="default");

    struct InputType
    {
      std::string data;
      bool isused;
    };
    typedef std::map<std::string, InputType > InputMap1d;
    typedef std::map<std::string, InputMap1d> InputMap2d;
    typedef std::map<std::string, InputMap2d> InputMap;
    InputMap inputlist;
    DataMap proflist;
    DataMap timelist;
    // std::vector<std::string> varnames;
    // std::vector<double> varvalues;
    std::string isused;
};
#endif
