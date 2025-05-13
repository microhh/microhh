import numpy as np

def get_data(case_dir):
    """
    Get `MOM/TKE/MASS` columns from `conservation.out`.
    """
    data_out = np.loadtxt(f'{case_dir}/conservation.out', skiprows=1)
    time = data_out[:,1]
    mom = data_out[:,7] / data_out[1,7]
    tke = data_out[:,8] / data_out[1,8]
    mass = data_out[:,9] / data_out[1,9]

    return time, mom, tke, mass


class Parse_conservation:
    def __init__(self, case_dir, experiment_name, scheme):
        """
        Parse momentum/tke/mass loss
        """

        dts = [1000, 500, 250, 125]
        self.N = len(dts)

        self.time = []
        self.mom = []
        self.tke = []
        self.mass = []

        self.mom_loss  = np.zeros(self.N)
        self.tke_loss  = np.zeros(self.N)
        self.mass_loss = np.zeros(self.N)

        for i,dt in enumerate(dts):
            time, mom, tke, mass = get_data(f'{case_dir}/{experiment_name}_{scheme}_dt{dt:04d}')

            self.time.append(time)
            self.mom.append(mom)
            self.tke.append(tke)
            self.mass.append(mass)

            self.mom_loss[i] = mom[-1] - mom[1]
            self.tke_loss[i] = tke[-1] - tke[1]
            self.mass_loss[i] = mass[-1] - mass[1]

        print(f'{scheme}: mom_loss = {self.mom_loss}')
        print(f'{scheme}: tke_loss = {self.tke_loss}')
        print(f'{scheme}: mass_loss = {self.mass_loss}')