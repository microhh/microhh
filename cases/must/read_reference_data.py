import numpy as np

class Read_2d_slice:
    def __init__(self, file_name):

        # Read data:
        f = np.loadtxt(file_name)
        self.x  = f[:,0 ]
        self.y  = f[:,1 ]
        self.z  = f[:,2 ]
        self.u  = f[:,3 ]
        self.v  = f[:,4 ]
        self.u2 = f[:,5 ]
        self.v2 = f[:,6 ]
        self.uv = f[:,7 ]
        self.u3 = f[:,8 ]
        self.v3 = f[:,9 ]
        self.u4 = f[:,10]
        self.v4 = f[:,11]

        # Unique x-and-y's:
        self.xu = np.unique(self.x)
        self.yu = np.unique(self.y)

class Read_tracer:
    def __init__(self, file_name):

        # Read data:
        f = np.loadtxt(file_name)
        self.i    = f[:,0]
        self.x    = f[:,1]
        self.y    = f[:,2]
        self.z    = f[:,3]
        self.c    = f[:,4]
        self.crms = f[:,5]
        self.lc   = f[:,6]
        self.c3   = f[:,7]
        self.c4   = f[:,8]

class Data_block:
    def __init__(self, block):
        self.x  = block[0,0 ]
        self.y  = block[0,1 ]
        self.z  = block[:,2 ]
        self.u  = block[:,3 ]
        self.w  = block[:,4 ]
        self.u2 = block[:,5 ]
        self.w2 = block[:,6 ]
        self.uw = block[:,7 ]
        self.u3 = block[:,8 ]
        self.w3 = block[:,9 ]
        self.u4 = block[:,10]
        self.w4 = block[:,10]

class Read_profile:
    def __init__(self, file_name):
        # Read data:
        f = np.loadtxt(file_name)

        x  = f[:,0 ]
        y  = f[:,1 ]

        # Get indices where profiles start/end
        ch = np.where((x[1:]-x[:-1] != 0) | (y[1:]-y[:-1] != 0))[0]
        ch = np.insert(ch, 0, -1)
        ch += 1
        ch = np.append(ch, x.size)

        # Save profiles in list
        self.data = []

        # Read all profiles
        for i in range(ch.size-1):
            self.data.append( Data_block(f[ch[i]:ch[i+1],:]) )

class Read_approach:
    def __init__(self, file_name):
        # Read data
        f = np.loadtxt(file_name)

        self.x  = f[:,0 ]
        self.y  = f[:,1 ]
        self.z  = f[:,2 ]

        self.u  = f[:,3 ]
        self.v  = f[:,4 ]
        self.w  = f[:,5 ]

        self.u2 = f[:,6 ]
        self.v2 = f[:,7 ]
        self.w2 = f[:,8 ]

        self.uv = f[:,9 ]
        self.uw = f[:,10]

        self.u3 = f[:,11]
        self.v3 = f[:,12]
        self.w3 = f[:,13]

        self.u4 = f[:,14]
        self.v4 = f[:,15]
        self.w4 = f[:,16]



if (__name__ == '__main__'):
    import matplotlib as mpl
    import matplotlib.pylab as pl
    from scipy.spatial import Voronoi, voronoi_plot_2d

    pl.close('all')

    z_ref = 7.29

    # -------------
    # Approach flow
    # -------------
    ap = Read_approach('reference_data/WT_approach_prof.txt')

    pl.figure()
    pl.subplot(221)
    pl.plot(ap.u, ap.z/z_ref, label='u')
    pl.plot(ap.v, ap.z/z_ref, label='v')
    pl.plot(ap.w, ap.z/z_ref, label='w')
    pl.legend()
    pl.ylabel('z/zref (-)')
    pl.xlabel('u_i/u_ref (-)')

    pl.subplot(222)
    pl.plot(ap.u2, ap.z/z_ref)
    pl.plot(ap.v2, ap.z/z_ref)
    pl.plot(ap.w2, ap.z/z_ref)
    pl.ylabel('z/zref (-)')
    pl.xlabel('u_i^`2/u_ref^2 (-)')

    pl.subplot(223)
    pl.plot(ap.u3, ap.z/z_ref)
    pl.plot(ap.v3, ap.z/z_ref)
    pl.plot(ap.w3, ap.z/z_ref)
    pl.ylabel('z/zref (-)')
    pl.xlabel('u_i^`3/u_ref^3 (-)')

    pl.subplot(224)
    pl.plot(ap.u4, ap.z/z_ref)
    pl.plot(ap.v4, ap.z/z_ref)
    pl.plot(ap.w4, ap.z/z_ref)
    pl.ylabel('z/zref (-)')
    pl.xlabel('u_i^`4/u_ref^4 (-)')

    pl.tight_layout()

    # -------------------------
    # Profiles zero degree case
    # -------------------------
    prf = Read_profile('reference_data/WT_00deg_prof.txt')

    pl.figure()
    for p in prf.data:
        pl.plot(p.u, p.z/z_ref)
    pl.xlabel('u/uref (-)')
    pl.ylabel('z/zref (-)')

    # -----------------
    # Tracer experiment
    # -----------------
    c = Read_tracer('reference_data/WT_45deg_tracer_1.28.txt')

    c.c[c.c<1e-4]=1e-4

    grid = np.vstack((c.x, c.y)).transpose()
    vor = Voronoi(grid)

    # find min/max values for normalization
    minima = np.min(np.log(c.c))
    maxima = np.max(np.log(c.c))

    # normalize chosen colormap
    norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
    mapper = pl.cm.ScalarMappable(norm=norm, cmap=pl.cm.jet)

    # plot Voronoi diagram, and fill finite regions with color mapped from speed value
    voronoi_plot_2d(vor, show_points=True, show_vertices=False, s=1, line_alpha=0.0)
    for r in range(len(vor.point_region)):
        region = vor.regions[vor.point_region[r]]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            pl.fill(*zip(*polygon), color=mapper.to_rgba(np.log(c.c[r])))

    pl.xlim(-95,95)
    pl.ylim(-100,10)
    pl.title('log(c)', loc='left')
    pl.xlabel('x (m)')
    pl.ylabel('y (m)')
