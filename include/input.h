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

#ifndef INPUT
#define INPUT

#include <map>
#include <string>
#include <vector>

// forward declaration to avoid circular dependency
class cmaster;

class cinput
{
  public:
    cinput(cmaster *);
    ~cinput();

    int readinput();
    int getItem(int *        , std::string, std::string, std::string);
    int getItem(int *        , std::string, std::string, std::string, int);
    int getItem(double *     , std::string, std::string, std::string);
    int getItem(double *     , std::string, std::string, std::string, double);
    int getItem(bool *       , std::string, std::string, std::string);
    int getItem(bool *       , std::string, std::string, std::string, bool);
    int getItem(std::string *, std::string, std::string, std::string);
    int getItem(std::string *, std::string, std::string, std::string, std::string);

    // list retrieval
    int getList(std::vector<int> *        , std::string, std::string, std::string);
    int getList(std::vector<std::string> *, std::string, std::string, std::string);

    int getProf(double *     , std::string, int size);
    int clear();
    int printUnused();

  private:
    cmaster *master;

    int readinifile ();
    int readproffile();

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
    int checkList(std::vector<std::string> *, std::string, std::string, std::string el="default");

    struct inputtype
    {
      std::string data;
      bool isused;
    };
    typedef std::map<std::string, inputtype> inputmap1d;
    typedef std::map<std::string, inputmap1d> inputmap2d;
    typedef std::map<std::string, inputmap2d> inputmap;
    inputmap inputlist;
    typedef std::map<std::string, std::vector<double> > profmap;
    profmap proflist;
    // std::vector<std::string> varnames;
    // std::vector<double> varvalues;
    std::string isused;
};
#endif
