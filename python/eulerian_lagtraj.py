from datetime import datetime
from lagtraj import DEFAULT_ROOT_DATA_PATH
from lagtraj.domain import download as domain_download
from lagtraj.forcings import create as forcing_create
from lagtraj.forcings import load as forcing_load
from lagtraj.trajectory import create as trajectory_create
import lagtraj
import os
import math
import numpy as np

import sys
sys.path.append('/usr/local/lib/python3.8/dist-packages/metpy')
sys.path.append('/home/girish/microhh_develop/microhh/python')
import netCDF4 as nc
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import math
import os
import microhh_tools as mht
import sched, time

float_type = 'f8'


def create_domain(domain_name,lat_min,lat_max,lon_min,lon_max): 
    print(f'Default lagtraj root path is: {DEFAULT_ROOT_DATA_PATH}')
    filen= str(DEFAULT_ROOT_DATA_PATH) + '/domains/' + domain_name + ".yaml"
    #f = open(filen,"a")
    f = open(filen,"w+")
    f.write(f"source    : era5\nversion   : 1.0.0\nlat_min   : {lat_min}\nlat_max   : {lat_max}\nlon_min   : {lon_min}\nlon_max   : {lon_max}\nlat_samp   : 0.1\nlon_samp   : 0.1")
    f.close()
    
def download_domain(domain_name,mm1,dd1,yyyy1,mm2,dd2,yyyy2):
    start_date=datetime(yyyy1,mm1,dd1)
    end_date=datetime(yyyy2,mm2,dd2)
    domain_download.download_named_domain(data_path=DEFAULT_ROOT_DATA_PATH,name=domain_name,start_date=start_date,end_date=end_date)
    args=[domain_name,f'{start_date}',f'{end_date}',"--data-path",f'{DEFAULT_ROOT_DATA_PATH}',"--retry-rate",'1']
    domain_download.cli(args=args)

def create_trajectory(domain_name,trajectory_name,lat,lon,mm,dd,yyyy,hrs_before,hrs_after,lagrangian): #start datetime "YYYY-MM-DDThh:ss"
    
    startdate_time=f"{yyyy}-{format(int(mm),'02d')}-{format(int(dd),'02d')}T00:00"
    print(startdate_time)
    filen= str(DEFAULT_ROOT_DATA_PATH) + '/trajectories/' + trajectory_name + ".yaml"


    f = open(filen,"w+")
    if lagrangian:
        f.write(f"trajectory_type : lagrangian\nvelocity_method : single_height_level\nvelocity_method_height : 500.\ndomain          : {domain_name}\n\
            version         : 1.0.1\nlat_origin      : {lat}\nlon_origin      : {lon}\ndatetime_origin : {startdate_time}\n\
backward_duration : PT{hrs_before}H\nforward_duration : PT{hrs_after}H\ntimestep        : domain_data")
    else:
        f.write(f"trajectory_type : eulerian\ndomain          : {domain_name}\nversion         : 1.0.1\nlat_origin      : {lat}\nlon_origin      : {lon}\ndatetime_origin : {startdate_time}\n\
backward_duration : PT{hrs_before}H\nforward_duration : PT{hrs_after}H\ntimestep        : domain_data")
    f.close()


    if os.path.exists(str(DEFAULT_ROOT_DATA_PATH)+'/trajectories/'+trajectory_name+'.nc'):
        os.remove(str(DEFAULT_ROOT_DATA_PATH)+'/trajectories/'+trajectory_name+'.nc')
    trajectory_create.main(data_path=DEFAULT_ROOT_DATA_PATH,trajectory_name=trajectory_name)


def create_forcing(domain_name,trajectory_name,conversion_type, runtime, surface_type, averaging_width):
    if surface_type=='ocean':
        gradient_method='regression'
        #averaging_width=2.0
    elif surface_type=='land':
        gradient_method='boundary'
        #averaging_width=1.0
    filen= str(DEFAULT_ROOT_DATA_PATH) + '/forcings/' + domain_name+ ".yaml"
    #f = open(filen,"a")
    f = open(filen,"w+")
    f.write(f"trajectory      : {trajectory_name}\nversion         : 1.0.1\ndomain          : {domain_name}\ngradient_method : {gradient_method}\n\
advection_velocity_sampling_method : domain_mean\nsampling_mask   : {surface_type}_only\naveraging_width : {averaging_width}\nlevels_method   : exponential\n\
levels_number   : 250\nlevels_dzmin    : 20.0\nlevels_ztop     : 72000.0")
    f.close()
    
    filen= str(DEFAULT_ROOT_DATA_PATH) + '/forcings/' + trajectory_name+ ".kpt.yaml"
    #f = open(filen,"a")
    f = open(filen,"w+")
    f.write(f"levels_method   : copy\nversion         : 1.0.0\nexport_format   : kpt\ncomment         : Forcing and initial conditions for Lagrangian case\n\
campaign        : EUREC4A\nsource_domain   : n/a\nreference       : n/a\nauthor          : s.j. boeing\nmodifications   : First version\ncase            : n/a\n\
adv_temp        : 1\nadv_theta       : 1\nadv_thetal      : 1\nadv_qv          : 1\nadv_qt          : 1\nadv_rv          : 1\nadv_rt          : 1\n\
rad_temp        : 0\nrad_theta       : 0\nrad_thetal      : 0\nforc_omega      : 0\nforc_w          : 1\nforc_geo        : 1\nnudging_u       : 0\nnudging_v       : 0\n\
nudging_temp    : 0\nnudging_theta   : 0\nnudging_thetal  : 0\nnudging_qv      : 0\nnudging_qt      : 0\nnudging_rv      : 0\nnudging_rt      : 0\nsurfaceType     : {surface_type}\n\
surfaceForcing  : ts\nsurfaceForcingWind               : z0_traj")
    f.close()
    if os.path.exists(str(DEFAULT_ROOT_DATA_PATH)+'/forcings/'+domain_name+'.nc'):
        os.remove(str(DEFAULT_ROOT_DATA_PATH)+'/forcings/'+domain_name+'.nc')
    forcing_defn=forcing_load.load_definition(DEFAULT_ROOT_DATA_PATH,forcing_name=domain_name)
    print(forcing_defn.name)
    if os.path.exists(str(DEFAULT_ROOT_DATA_PATH)+'/forcings/'+domain_name+'.kpt.nc'):
        os.remove(str(DEFAULT_ROOT_DATA_PATH)+'/forcings/'+domain_name+'.kpt.nc')
    forcing_create.main(data_path=DEFAULT_ROOT_DATA_PATH,forcing_defn=forcing_defn,conversion_name=conversion_type)
    create_microhhforcing(domain_name,tstart=runtime)

def create_microhhforcing(domain_name,tstart):
    basename = str(DEFAULT_ROOT_DATA_PATH)+'/forcings/'+domain_name
    netcdf_path = basename+".kpt.nc"
    evergreen=0.7;
    all_data = xr.open_dataset(netcdf_path,decode_times=False)
    time = all_data['time'].values[0:];
    select_arr=(time>=time[tstart])
    qadv_un=all_data['qadv'].values[select_arr,:];             qadv_un=np.flip(qadv_un, axis=1);
    tadv_un=all_data['tadv'].values[select_arr,:];             tadv_un=np.flip(tadv_un, axis=1);
    uadv_un=all_data['uadv'].values[select_arr,:];             uadv_un=np.flip(uadv_un, axis=1);
    vadv_un=all_data['vadv'].values[select_arr,:];             vadv_un=np.flip(vadv_un, axis=1);
    time = all_data['time'].values[select_arr];
    qt_un=all_data['q'].values[select_arr,:];                  qt_un=np.flip(qt_un, axis=1);
    ql_un=all_data['ql'].values[select_arr,:];                 ql_un=np.flip(ql_un, axis=1);
    u_un=all_data['u'].values[select_arr,:];                   u_un=np.flip(u_un, axis=1);
    v_un=all_data['v'].values[select_arr,:];                   v_un=np.flip(v_un, axis=1);
    ug_un=all_data['ug'].values[select_arr,:];                 ug_un=np.flip(ug_un, axis=1);
    vg_un=all_data['vg'].values[select_arr,:];                 vg_un=np.flip(vg_un, axis=1);
    zun=all_data['zf'].values[select_arr,:];             zun=np.flip(zun, axis=1);
    pres_un=all_data['pres'].values[select_arr,:];       pres_un=np.flip(pres_un, axis=1);
    T_un=all_data['t'].values[select_arr,:];                   T_un=np.flip(T_un, axis=1);
    pres0=all_data['ps'].values[select_arr];
    sst=all_data['t_skin'].values[select_arr];
    qs=all_data['q_skin'].values[select_arr];
    z0m=all_data['mom_rough'].values[select_arr];
    z0h=all_data['heat_rough'].values[select_arr];
    H=all_data['sfc_sens_flx'].values[select_arr];
    albedo=all_data['albedo'].values[select_arr];
    LE = all_data['sfc_lat_flx'].values[select_arr];
    omega_un = all_data['omega'].values[select_arr,:];         omega_un=np.flip(omega_un, axis=1);
    o3_un = all_data['o3'].values[select_arr,:];               o3_un=np.flip(o3_un, axis=1);
    time = all_data['time'].values[select_arr];
    mean_height=all_data['orog'].values[:]; 
    lon = all_data['lon'].values[select_arr]; 
    lat = all_data['lat'].values[select_arr];
    lat = all_data['lat'].values[select_arr];
    h_soil = all_data['h_soil'].values[:];
    q_soil = all_data['q_soil'].values[select_arr,:];
    t_soil = all_data['t_soil'].values[select_arr,:];
    low_veg_lai = all_data['low_veg_lai'].values[select_arr];
    high_veg_lai = all_data['high_veg_lai'].values[select_arr];
    low_veg_cover = all_data['low_veg_cover'].values[select_arr];
    high_veg_cover = all_data['high_veg_cover'].values[select_arr];

    # set the height
    z_new=np.zeros(300)
    dz=10
    z_abs=10;
    for i in range(0,z_new.size): 
            if i<12:
                dz=dz+int(round(0.1*dz,0));
            elif i==14:
                z_abs=300;
                dz=40;
            elif z_abs>5000:
                dz=dz+int(round(0.1*dz,0));
            z_new[i] = z_abs;
            z_abs = z_abs+dz
    z_end_ind=np.nonzero((z_new>29000))[0][0]
    z=z_new[0:z_end_ind+1]
    kmax=z.size
    ############################## Declare input variables to nc input and constants ##################################

    sat_r = np.zeros(time.size)
    qt_bot = np.zeros(time.size)
    sbotthl = np.zeros(time.size)

    u = np.zeros((time.size, kmax))
    v = np.zeros(np.shape(u))
    ugeo = np.zeros(np.shape(u))
    vgeo = np.zeros(np.shape(u))
    qt = np.zeros(np.shape(u))
    ql = np.zeros(np.shape(u))
    qadv = np.zeros(np.shape(u))
    tadv = np.zeros(np.shape(u))
    uadv = np.zeros(np.shape(u))
    vadv = np.zeros(np.shape(u))
    th = np.zeros(np.shape(u))
    thl   = np.zeros(np.shape(u))
    thlls = np.zeros(np.shape(u))
    qtls = np.zeros(np.shape(u))
    w   = np.zeros(np.shape(u))
    uls = np.zeros(np.shape(u))
    vls = np.zeros(np.shape(u))
    wls = np.zeros((time.size, kmax+1))
    pres = np.zeros(np.shape(u))
    omega = np.zeros(np.shape(u))
    o3_f = np.zeros(np.shape(u))
    T = np.zeros(np.shape(u))
    nudge_factor = np.zeros(np.shape(u))
    th_diff = np.zeros(time.size)
    qt_diff = np.zeros(time.size)
    U = np.zeros(time.size)

    cp  = 1005.
    Lv  = 2.5e6
    Rd  = 287.
    tau = 10800;

    ######################## Radiation Calculation and NC input ##################################

    nc_file = nc.Dataset(basename+"_input.nc", mode="w", datamodel="NETCDF4", clobber=True)
    
    z_top = 70.e3
    dz = 500.
    z_rad  = np.arange(dz/2, z_top, dz)
    zh_rad = np.arange(   0, z_top-dz/2, dz)
    zh_rad = np.append(zh_rad, z_top)
    print(mean_height, zun.shape)
    #print(zun[0,:])
    if np.isnan(mean_height[0]):
        mean_height[0]=0.
    zun_rad=zun[0,:]-mean_height[0]
    interp_rad=(np.logical_not(np.isnan(zun_rad)))
    p_lay = np.interp(z_rad,zun_rad[interp_rad],pres_un[0,interp_rad])
    p_lev = np.interp(zh_rad,zun_rad[interp_rad],pres_un[0,interp_rad])
    T_lay = np.interp(z_rad,zun_rad[interp_rad],T_un[0,interp_rad])
    T_lev = np.interp(zh_rad,zun_rad[interp_rad],T_un[0,interp_rad])

    co2 =  348.e-6
    ch4 = 1650.e-9
    n2o =  306.e-9
    n2 = 0.7808
    o2 = 0.2095

    nc_group_rad = nc_file.createGroup("radiation")

    nc_group_rad.createDimension("lay", p_lay.size)
    nc_group_rad.createDimension("lev", p_lev.size)

    nc_z_lay = nc_group_rad.createVariable("z_lay", float_type, ("lay"))
    nc_z_lev = nc_group_rad.createVariable("z_lev", float_type, ("lev"))
    nc_z_lay[:] = z_rad [:]
    nc_z_lev[:] = zh_rad[:]

    nc_p_lay = nc_group_rad.createVariable("p_lay", float_type, ("lay"))
    nc_p_lev = nc_group_rad.createVariable("p_lev", float_type, ("lev"))
    nc_p_lay[:] = p_lay[:]
    nc_p_lev[:] = p_lev[:]

    nc_T_lay = nc_group_rad.createVariable("t_lay", float_type, ("lay"))
    nc_T_lev = nc_group_rad.createVariable("t_lev", float_type, ("lev"))
    nc_T_lay[:] = T_lay[:]
    nc_T_lev[:] = T_lev[:]

    nc_CO2 = nc_group_rad.createVariable("co2", float_type)
    nc_CH4 = nc_group_rad.createVariable("ch4", float_type)
    nc_N2O = nc_group_rad.createVariable("n2o", float_type)
    nc_O3  = nc_group_rad.createVariable("o3" , float_type, ("lay"))
    nc_H2O = nc_group_rad.createVariable("h2o", float_type, ("lay"))
    nc_N2  = nc_group_rad.createVariable("n2" , float_type)
    nc_O2  = nc_group_rad.createVariable("o2" , float_type)

    nc_CFC11 = nc_group_rad.createVariable("cfc11", float_type)
    nc_CFC12 = nc_group_rad.createVariable("cfc12", float_type)
    nc_CFC22 = nc_group_rad.createVariable("cfc22", float_type)
    nc_CCL4  = nc_group_rad.createVariable("ccl4" , float_type)

    nc_CO2[:] = co2
    nc_CH4[:] = ch4
    nc_N2O[:] = n2o
    #nc_O3 [:] = o3 [:]
    #nc_H2O[:] = h2o[:]
    nc_N2 [:] = n2
    nc_O2 [:] = o2

    nc_CFC11[:] = 0.
    nc_CFC12[:] = 0.
    nc_CFC22[:] = 0.
    nc_CCL4 [:] = 0.

    qt_rad=np.zeros(z_rad.size); o3_rad=np.zeros(z_rad.size);
    qt_rad[:] = np.interp(z_rad,zun_rad[interp_rad],qt_un[0,interp_rad])
    o3_rad[:] = np.interp(z_rad,zun_rad[interp_rad],o3_un[0,interp_rad])

    xm_air = 28.97; xm_h2o = 18.01528
    h2o=qt_rad*xm_air/xm_h2o
    nc_H2O[:] = h2o[:]
    nc_O3[:] = o3_rad[:]
    ######################## Calculation of variables ############################################

    zh = 0.5*(z[:-1] + z[1:])
    zh = np.append(0., zh)
    zh = np.append(zh, z.size)

    time = time - time[0];


    for n in range(0,time.size):
        if np.isnan(mean_height[n]):
            mean_height[n]=0.
        zun[n,:]=zun[n,:]-mean_height[n]
        interp_arr=(np.logical_not(np.isnan(zun[n,:])))
        qt[n,:] = np.interp(z,zun[n,interp_arr],qt_un[n,interp_arr])
        ql[n,:] = np.interp(z,zun[n,interp_arr],ql_un[n,interp_arr])
        u[n,:] = np.interp(z,zun[n,interp_arr],u_un[n,interp_arr])
        v[n,:] = np.interp(z,zun[n,interp_arr],v_un[n,interp_arr])
        ugeo[n,:] = np.interp(z,zun[0,interp_arr],ug_un[0,interp_arr])
        vgeo[n,:] = np.interp(z,zun[0,interp_arr],vg_un[0,interp_arr])
        omega[n,:] = np.interp(z,zun[n,interp_arr],omega_un[n,interp_arr])
        o3_f[n,:] = np.interp(z,zun[n,interp_arr],o3_un[n,interp_arr])
        pres[n,:] = np.interp(z,zun[n,interp_arr],pres_un[n,interp_arr])
        T[n,:] = np.interp(z,zun[n,interp_arr],T_un[n,interp_arr])
        qadv[n,:] = np.interp(z,zun[n,interp_arr],qadv_un[n,interp_arr])
        tadv[n,:] = np.interp(z,zun[n,interp_arr],tadv_un[n,interp_arr])
        uadv[n,:] = np.interp(z,zun[n,interp_arr],uadv_un[n,interp_arr])
        vadv[n,:] = np.interp(z,zun[n,interp_arr],vadv_un[n,interp_arr])
        qt[n,:] = np.interp(z,zun[n,interp_arr],qt_un[n,interp_arr])
        u[n,:] = np.interp(z,zun[n,interp_arr],u_un[n,interp_arr])
        v[n,:] = np.interp(z,zun[n,interp_arr],v_un[n,interp_arr])
        ugeo[n,:] = np.interp(z,zun[n,interp_arr],ug_un[n,interp_arr])
        vgeo[n,:] = np.interp(z,zun[n,interp_arr],vg_un[n,interp_arr])

    ug = u; vg = v;
    p_sbot = pres[:,0];
    nudge_factor[:,:]=1./tau

    for n in range(0,time.size):
        sat_r = mpcalc.saturation_mixing_ratio(p_sbot[n] * units.pascal , sst[n]* units.kelvin)
        qt_bot[n] = 0.981 * mpcalc.specific_humidity_from_mixing_ratio(sat_r)
        
        for k in range(0,kmax):
            w[n,k] = mpcalc.vertical_velocity(omega[n,k] * units.pascal / units.second, pres[n,k] * units.pascal, T[n,k] * units.kelvin) / (units.meter / units.second)
            th[n,k] = mpcalc.potential_temperature(pres[n,k] * units.pascal, T[n,k] * units.kelvin) / units.kelvin
            thl[n,k] = th[n,k] - (th[n,k]/T[n,k]) * (Lv/cp) * (ql[n,k]/(1-qt[n,k]))
    fc_cal = mpcalc.coriolis_parameter(np.mean(lat)*units.degrees) * units.second
    for n in range(0,time.size):
        wls[n,:] = np.interp(zh,z,w[n,:])

    ### Fluxes ###
    rhosurf = p_sbot / (Rd * thl[:,0] * (1. + 0.61 * qt[:,0]))
    lh_flx = -LE / (rhosurf * Lv) #J/m2s / (J/m3) --> m/s
    sh_flx = -H / (rhosurf * cp) # K m/s
    
    ######################################### Land Surface Model #######################################
    def link(f1, f2):
        """
        Create symbolic link from `f1` to `f2`, if `f2` does not yet exist.
        """
        if os.path.islink(f2):
            os.remove(f2)
        if os.path.exists(f1):
            os.symlink(f1, f2)
        else:
            raise Exception('Source file {} does not exist!'.format(f1))

    def add_nc_var(name, dims, nc, data):
        """
        Create NetCDF variables and set values.
        """
        if dims is None:
            var = nc.createVariable(name, np.float64)
        else:
            var = nc.createVariable(name, np.float64, dims)
        var[:] = data

    ############################## write the data to a file ############################################

    nc_file.createDimension("z", kmax)
    nc_z = nc_file.createVariable("z", float_type, ("z"))
    nc_z    [:] = z    [:]

    nc_file.createDimension("zh", kmax+1)
    nc_zh = nc_file.createVariable("zh", float_type, ("zh"))
    nc_zh[:] = zh[:]

    nc_file.createDimension("time_ls", time.size)
    nc_time_ls = nc_file.createVariable("time_ls", float_type, ("time_ls"))
    nc_time_ls [:] = time [:]

    nc_file.createDimension("time_surface", time.size)
    nc_time_surface = nc_file.createVariable("time_surface", float_type, ("time_surface"))
    nc_time_surface [:] = time [:]

    nc_group_init = nc_file.createGroup("init");
    nc_group_timedep = nc_file.createGroup("timedep")

    ##### initial conditions ############
    nc_thl   = nc_group_init.createVariable("thl"   , float_type, ("z"))
    nc_qt    = nc_group_init.createVariable("qt"    , float_type, ("z"))
    nc_u     = nc_group_init.createVariable("u"     , float_type, ("z"))
    nc_ugeo  = nc_group_init.createVariable("u_geo" , float_type, ("z"))
    nc_v     = nc_group_init.createVariable("v"     , float_type, ("z"))
    nc_vgeo  = nc_group_init.createVariable("v_geo" , float_type, ("z"))
    nc_wls  = nc_group_init.createVariable("w_ls" , float_type, ("zh"))
    nc_qtls  = nc_group_init.createVariable("qt_ls" , float_type, ("z"))
    nc_thlls  = nc_group_init.createVariable("thl_ls" , float_type, ("z"))

    nc_CO2 = nc_group_init.createVariable("co2", float_type)
    nc_CH4 = nc_group_init.createVariable("ch4", float_type)
    nc_N2O = nc_group_init.createVariable("n2o", float_type)
    nc_O3  = nc_group_init.createVariable("o3" , float_type, ("z"))
    nc_H2O = nc_group_init.createVariable("h2o", float_type, ("z"))
    nc_N2  = nc_group_init.createVariable("n2" , float_type)
    nc_O2  = nc_group_init.createVariable("o2" , float_type)

    nc_CFC11 = nc_group_init.createVariable("cfc11", float_type)
    nc_CFC12 = nc_group_init.createVariable("cfc12", float_type)
    nc_CFC22 = nc_group_init.createVariable("cfc22", float_type)
    nc_CCL4  = nc_group_init.createVariable("ccl4" , float_type)


    ###### forcing conditions ############
    nc_group_timedep.createDimension("time_ls", time.size)
    nc_time_ls = nc_group_timedep.createVariable("time_ls", float_type, ("time_ls"))
    nc_time_ls [:] = time [:]

    nc_group_timedep.createDimension("time_nudge", time[:].size)
    nc_time_nudge = nc_group_timedep.createVariable("time_nudge", float_type, ("time_nudge"))
    nc_time_nudge [:] = time [:]

    nc_group_timedep.createDimension("time_surface", time.size)
    nc_time_surface = nc_group_timedep.createVariable("time_surface", float_type, ("time_surface"))
    nc_time_surface [:] = time [:]

    nc_u_ls   = nc_group_timedep.createVariable("u_ls" , float_type, ("time_ls","z"))
    nc_v_ls   = nc_group_timedep.createVariable("v_ls" , float_type, ("time_ls","z"))
    nc_u_g = nc_group_timedep.createVariable("u_geo", float_type, ("time_ls", "z"))
    nc_v_g = nc_group_timedep.createVariable("v_geo", float_type, ("time_ls", "z"))
    nc_w_ls   = nc_group_timedep.createVariable("w_ls" , float_type, ("time_ls","zh"))
    nc_thl_ls = nc_group_timedep.createVariable("thl_ls" , float_type, ("time_ls","z"))
    nc_qt_ls  = nc_group_timedep.createVariable("qt_ls" , float_type, ("time_ls","z")) 


    ###### nudge conditions ##############
    
    nc_nudge_factor = nc_group_init.createVariable("nudgefac", float_type, ("z"))
    nc_u_nudge = nc_group_timedep.createVariable(
        "u_nudge", float_type, ("time_nudge", "z"))
    nc_v_nudge = nc_group_timedep.createVariable(
        "v_nudge", float_type, ("time_nudge", "z"))
    nc_thl_nudge = nc_group_timedep.createVariable(
        "thl_nudge", float_type, ("time_nudge", "z"))
    nc_qt_nudge = nc_group_timedep.createVariable(
        "qt_nudge", float_type, ("time_nudge", "z"))
    ###### time dependent bottom conditions ####### 
    nc_thl_sbot = nc_group_timedep.createVariable("thl_sbot", float_type, ("time_surface"))
    nc_qt_sbot = nc_group_timedep.createVariable("qt_sbot", float_type, ("time_surface"))
    nc_p_sbot = nc_group_timedep.createVariable("p_sbot", float_type, ("time_surface"))



    nc_thl  [:] = thl  [0,:]
    nc_qt   [:] = qt   [0,:]
    nc_u    [:] = u    [0,:]
    nc_ugeo [:] = ug   [0,:]
    nc_v    [:] = v    [0,:]
    nc_vgeo [:] = vg   [0,:]
    nc_wls  [:] = wls  [0,:]
    nc_qtls [:] = qadv [0,:]
    nc_thlls[:] = tadv [0,:]

    nc_u_g  [:, :] = ug  [:, :]
    nc_v_g  [:, :] = vg  [:, :]
    nc_u_ls  [:, :] = uadv  [:, :]
    nc_v_ls  [:, :] = vadv  [:, :]
    nc_w_ls  [:, :] = wls  [:, :]
    nc_thl_ls[:, :] = tadv [:, :]
    nc_qt_ls [:, :] = qadv [:, :]

    nc_CO2[:] = co2
    nc_CH4[:] = ch4
    nc_N2O[:] = n2o
    nc_O3 [:] = o3_f[0,:]
    nc_H2O[:] = qt[0,:] * xm_air/xm_h2o
    nc_N2 [:] = n2
    nc_O2 [:] = o2

    nc_CFC11[:] = 0.
    nc_CFC12[:] = 0.
    nc_CFC22[:] = 0.
    nc_CCL4 [:] = 0.

    nc_thl_sbot[:] = sh_flx[:]
    nc_qt_sbot[:] = lh_flx[:]
    nc_p_sbot[:] = p_sbot[:]

    #### if nudge #####
    nc_u_nudge[:, :] = u[:, :]
    nc_v_nudge[:, :] = v[:, :]
    nc_thl_nudge[:, :] = thl[:, :]
    nc_qt_nudge[:, :] = qt[:, :]
    nc_nudge_factor[:] = nudge_factor[0, :]


    

    ################################## update ini file ########################
    ini = mht.Read_namelist(basename+'.ini.base')

    ini['grid']['ktot'] = kmax
    ini['grid']['zsize'] = z_new[kmax]
    ini['thermo']['pbot'] = p_sbot[0]
    ini['radiation']['lat'] = np.mean(lat)
    ini['radiation']['lon'] = np.mean(lon)
    ini['radiation']['sfc_alb_dir'] = np.mean(albedo)
    ini['radiation']['sfc_alb_dif'] = np.mean(albedo)
    ini['force']['fc'] = fc_cal.magnitude
    ini['boundary']['z0m'] = z0m[0]
    ini['boundary']['z0h'] = z0h[0]
    
    if ini['boundary']['swboundary'] == 'surface_lsm':
        use_htessel = True
    else:
        use_htessel = False
    
    
    if use_htessel:
        
        type_soil=3;
        root_frac = np.zeros(np.shape(h_soil))
        #root_frac = [0.244760790777786, 0.409283067913477, 0.307407403941806, 0.0385487373669315]

        if domain_name=='SEUS':
            a=6.706;
            b=2.175;
        elif domain_name=='SGP':
            a=5.558;
            b=2.614;
        root_frac=1-0.5*(np.exp(-a*h_soil)+np.exp(-b*h_soil))
        for n in range(len(h_soil)-1,0,-1):
            root_frac[n]=root_frac[n]-(root_frac[n-1])

        # link('/home/girish/microhh_develop/microhh/misc/van_genuchten_parameters.nc', 'van_genuchten_parameters.nc')
        nc_group_soil = nc_file.createGroup("soil")
        nc_group_soil.createDimension('z', h_soil.size)
        index_soil = np.ones_like(h_soil)*int(type_soil-1)

        add_nc_var('z', ('z'), nc_group_soil, -h_soil[::-1])

        add_nc_var('theta_soil', ('z'), nc_group_soil, q_soil[0,::-1])
        add_nc_var('t_soil', ('z'), nc_group_soil, t_soil[0,::-1])
        add_nc_var('index_soil', ('z'), nc_group_soil, index_soil)
        add_nc_var('root_frac', ('z'), nc_group_soil, root_frac[::-1])

          
        ini['boundary']['sbcbot'] = 'dirichlet'
        ini['land_surface']['swhomogeneous'] = True
        ini['boundary']['swconstantz0'] = True
        
        if domain_name=='SEUS':
            gD_hv=evergreen*0.0003+(1-evergreen)*0.0013;
            #rs_highveg=evergreen*500+(1-evergreen)*240;
            lai=high_veg_cover[0]*high_veg_lai[0]+low_veg_cover[0]*low_veg_lai[0]
            ini['land_surface']['lai'] = 5.5
            ini['land_surface']['gD'] = high_veg_cover[0]*gD_hv
            ini['land_surface']['lambda_stable'] = high_veg_cover[0]*20+low_veg_cover[0]*10
            ini['land_surface']['lambda_unstable'] = high_veg_cover[0]*15+low_veg_cover[0]*10
            ini['land_surface']['c_veg'] = 0.99
            ini['radiation']['emis_sfc'] = 0.97
            ini['land_surface']['rs_veg_min'] = 180
        elif domain_name=='SGP':
            lai=high_veg_cover[0]*high_veg_lai[0]+low_veg_cover[0]*low_veg_lai[0]
            ini['land_surface']['lai'] = lai
            ini['land_surface']['gD'] = 0
            ini['land_surface']['lambda_stable'] = 10
            ini['land_surface']['lambda_unstable'] = 10
            ini['land_surface']['c_veg'] = high_veg_cover[0]*0.6+low_veg_cover[0]*1

    else:
        ini['boundary']['swboundary'] = 'surface'
        ini['boundary']['sbcbot'] = 'flux'
        ini['boundary']['swtimedep'] = True
        ini['boundary']['timedeplist'] = ['thl_sbot','qt_sbot']
        ini['boundary']['sbot[thl]'] = sh_flx[0]
        ini['boundary']['sbot[qt]'] = lh_flx[0]
        ini['boundary']['stop[qt]'] = 0
        ini['boundary']['stop[thl]'] = 0
        
    ini.save(basename+'.ini', allow_overwrite=True)
    nc_file.close()



# domain = input('Which domain is the forcing needed for? (SGP, SEUS, others): \n')
# if domain!='SGP' and domain!='SEUS' and domain!='others':
#     domain = input('Please enter one from the options --> (SGP, SEUS, others): \n')
# if domain=='others':
#     domain_name = input('Please provide the domain name: \n')
#     create_new_domain = input('Do you want to create a new domain? (y,n):')
#     if create_new_domain=='y':
#         latd,lond = np.float64(input('Please enter the Latitude and Longitude of the domain: \n').split())
#         create_domain(f'{domain_name}',lat_min=math.floor(latd-1),lat_max=math.ceil(latd),lon_min=math.floor(lond-1),lon_max=math.ceil(lond))
# elif domain=='SGP' or domain=='SEUS':
#     domain_name=domain

# domaindown = input('Do you want to download domain data? (y,n):')
# if domaindown=='y':
#     start_date=input('Please enter the start date to download domain (YYYY,MM,DD): \n').split(",")
#     end_date=input('Please enter the end date to download domain (YYYY,MM,DD): \n').split(",")
#     download_domain(domain_name,mm1=int(start_date[1]),dd1=int(start_date[2]),yyyy1=int(start_date[0]),mm2=int(end_date[1]),dd2=int(end_date[2]),yyyy2=int(end_date[0]))   

# forcingcreate = input('Do you want to create forcing? (y,n):')

domain_name = 'BOMEX'
latd = 11.5
lond = -55.5
start_date = ['1969','06','22']
end_date = ['1969','06','23']
hours_before = 0
hours_after = 24
surface_type='ocean'
averaging_width=2.0




islagrangian=False
create_new_domain = False
download_domain_data = False
create_new_forcing = True

domain_name='SEUS'
start_date = ['2015','08','09']
end_date = ['2015','08','13']
hours_before = 0
hours_after = 120
averaging_width=1.0

if create_new_domain:
    create_domain(f'{domain_name}',lat_min=math.floor(latd-1),lat_max=math.ceil(latd),lon_min=math.floor(lond-1),lon_max=math.ceil(lond))

if download_domain_data:
    
    download_domain(domain_name,mm1=int(start_date[1]),dd1=int(start_date[2]),yyyy1=int(start_date[0]),\
        mm2=int(end_date[1]),dd2=int(end_date[2]),yyyy2=int(end_date[0]))   


if create_new_forcing:
    # start_date = input('Please enter the start date of the trajectory (YYYY,MM,DD): \n').split(",") 
    # hours_after = int(input('Please enter the hours after start date till which trajctory is needed for: \n'))

    
    if domain_name=='SGP':
        create_trajectory(domain_name,domain_name,lat=36.5,lon=-97.5, \
        mm=int(start_date[1]),dd=int(start_date[2]),yyyy=int(start_date[0]),hrs_before=hours_before,hrs_after=hours_after, lagrangian=islagrangian)
        create_forcing(domain_name,domain_name,conversion_type="kpt", runtime = 0, surface_type='land', averaging_width=averaging_width)
    elif domain_name=='SEUS':
        create_trajectory(domain_name,domain_name,lat=34.25,lon=-87.25, \
        mm=int(start_date[1]),dd=int(start_date[2]),yyyy=int(start_date[0]),hrs_before=hours_before,hrs_after=hours_after, lagrangian=islagrangian)
        create_forcing(domain_name,domain_name,conversion_type="kpt", runtime = 0, surface_type='land', averaging_width=averaging_width)
    else:
        create_trajectory(domain_name,domain_name,lat=latd,lon=lond, \
        mm=int(start_date[1]),dd=int(start_date[2]),yyyy=int(start_date[0]),hrs_before=hours_before,hrs_after=hours_after, lagrangian=islagrangian)
        create_forcing(domain_name,domain_name,conversion_type="kpt", runtime = 0, surface_type=surface_type, averaging_width=averaging_width)
    