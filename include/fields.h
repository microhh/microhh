/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

#ifndef FIELDS
#define FIELDS
#include <map>
#include <vector>
#include "field3d.h"

// forward declarations to reduce compilation time
class Master;
class Input;
class Model;
class Grid;
class Stats;
struct mask;

typedef std::map<std::string, Field3d *> fieldmap;

class Fields
{
  public:
    // functions
    Fields(Model *, Input *);
    ~Fields();

    void init();
    void create(Input *);

    int exec();
    int getmask(Field3d *, Field3d *, mask *);
    int execstats(mask *);

    int initmomfld(Field3d*&, Field3d*&, std::string, std::string, std::string);
    int initpfld(std::string, std::string, std::string);
    int initdfld(std::string, std::string, std::string);
    
    void save(int);
    void load(int);

    double checkmom ();
    double checktke ();
    double checkmass();

    int setcalcprofs(bool);

    void execcross();

    // 3d fields for momentum
    Field3d *u;
    Field3d *v;
    Field3d *w;

    Field3d *ut;
    Field3d *vt;
    Field3d *wt;

    // maps of 3d fields
    fieldmap a;
    fieldmap ap;
    fieldmap at;

    fieldmap mp;
    fieldmap mt;

    fieldmap s;
    fieldmap sd;
    fieldmap sp;
    fieldmap st;

    // reference density
    double *rhoref;
    double *rhorefh;

    // TODO remove these to and bring them to diffusion model
    double visc;

    // GPU functions and variables
    enum OffsetType {Offset, NoOffset};

    int prepareDevice();
    int forwardDevice();
    void forward3DFieldDevice (double *, double *, OffsetType);
    void forward2DFieldDevice (double *, double *, OffsetType);
    void forward1DFieldDevice (double *, double *, int);
    void backward3DFieldDevice(double *, double *, OffsetType);
    void backward2DFieldDevice(double *, double *, OffsetType);
    void backward1DFieldDevice(double *, double *, int);
    int backwardDevice();
    int clearDevice();

    double *rhoref_g;
    double *rhorefh_g;
    
  private:
    // variables
    Model  *model;
    Grid   *grid;
    Master *master;
    Stats  *stats;

    bool calcprofs;

    // cross sections
    std::vector<std::string> crosslist;      // List with all crosses from ini file
    // Cross sections split per type
    std::vector<std::string> crosssimple;
    std::vector<std::string> crosslngrad;   
    std::vector<std::string> crossbot;
    std::vector<std::string> crosstop;
    std::vector<std::string> crossfluxbot;
    std::vector<std::string> crossfluxtop;
    int checkaddcross(std::string, std::string, std::vector<std::string> *, std::vector<std::string> *);

    // masks
    int calcmaskwplus(double *, double *, double *, int *, int *, int *, double *);
    int calcmaskwmin (double *, double *, double *, int *, int *, int *, double *);

    // perturbations
    double rndamp;
    double rndz;
    double rndexp;
    double vortexamp;
    int vortexnpair;
    std::string vortexaxis;
    
    // functions
    double calcmom_2nd(double *, double *, double *, double *);
    double calctke_2nd(double *, double *, double *, double *);
    int addmeanprofile(Input *, std::string, double *, double);
    int randomnize(Input *, std::string, double *);
    int addvortexpair(Input* inputin);
    double calcmass(double *, double *);

    // statistics
    double *umodel;
    double *vmodel;
};
#endif
