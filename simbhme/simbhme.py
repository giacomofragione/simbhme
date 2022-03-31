'''
simbhme: computing cosmological merger rates for supermassive black holes - intermediate-mass black holes binaries
https://github.com/giacomofragione/smbh_imbh
'''

import numpy as np
import os
import sys

import astropy.units as u
from astropy.cosmology import Planck15, z_at_value

import random
from datetime import datetime
random.seed(datetime.now())


########################################################3

__author__ = "Giacomo Fragione"
__license__ = "MIT"
__version__ = "1.0.1"
__email__ = "giacomo.fragione@northwestern.edu"
this_module='simbhme'

########################################################3

#Defaults values
defaults={ 'directory' : os.path.dirname(__file__),
           'mgal_min' : 10 ** 8.5,
           'mgal_max' : 10 ** 10.75,
           'mstar' : 10 ** 11.14,
           'astar' : 1.43,
           'galaxy_type' : 'late',
           'mnsc_min' : 1e5,
           'mnsc_max' : 1e8,
           'mgc_min' : 1e5,
           'mgc_max' : 1e7,
           'slope_gc' : 2,
           'fimbh' : 0.001,
           'fgcimbh' : 1.0,
           'tdelay' : 1,
           'scale' : 0,
           'z_gc' : 3.2,
           'sigma_gc' : 1.5}


class smbh_imbh_mergers(object):
    '''
    Compute cosmological merger rates for supermassive black holes - intermediate-mass black holes binaries

    Usage:
        p=sim.smbh_imbh_mergers()
        p(nsample)

    Parameters:
        nsample # number of galaxies to sample
        directory # directory
        mnsc_min # minimum NSC mass
        mnsc_max # maximum NSC mass
        mgc_min # minimum GC mass
        mgc_max # maximum GC mass
        slope_gc # slope GC distribution
        z_gc # mean redshift GC formation
        sigma_gc # sigma redshift GC formation
        mgal_min # minimum galaxy mass
        mgal_max # maximum galaxy mass late
        mstar  # parameter press-schechter function
        astar  # parameter press-schechter function
        galaxy_type # galaxy type
        fimbh # IMBH mass over cluster mass
        fgcimbh # fraction of clusters with IMBH
        tdelay # delay time in Gyr
        scale # =1 determine minimum GC mass with scaling from Fahrion+ 2021

    Returns:
        akk # number of galaxy
        amgal # galaxy mass
        amnsc # NSC mass
        arnsc # half-mass radius NSC
        amsmbh # SMBH mass
        arhsmbh # influence radius SMBH
        afout # NSC mass from accreted GCs
        amgc # GC mass
        azform # formation redshit GC
        atform # formation time GC
        amimbh # IMBH mass
        attot # total time for SMBH-IMBH merger
        aii # cluster index
        akkk # number of runs to meet criteria of GC mass generation
    '''


    def __init__(self,  directory=defaults['directory'],
                        mgal_min=defaults['mgal_min'],
                        mgal_max=defaults['mgal_max'],
                        mstar=defaults['mstar'],
                        astar=defaults['astar'],
                        galaxy_type=defaults['galaxy_type'],
                        mnsc_min=defaults['mnsc_min'],
                        mnsc_max=defaults['mnsc_max'],
                        mgc_min=defaults['mgc_min'],
                        mgc_max=defaults['mgc_max'],
                        slope_gc=defaults['slope_gc'],
                        fimbh=defaults['fimbh'],
                        fgcimbh=defaults['fgcimbh'],
                        tdelay=defaults['tdelay'],
                        scale=defaults['scale'],
                        z_gc=defaults['z_gc'],
                        sigma_gc=defaults['sigma_gc']):

        self.directory = directory
        self.mgal_min = mgal_min           
        self.mgal_max = mgal_max
        self.mstar = mstar
        self.astar = astar
        self.galaxy_type = galaxy_type
        self.mnsc_min = mnsc_min
        self.mnsc_max = mnsc_max
        self.mgc_min = mgc_min
        self.mgc_max = mgc_max
        self.slope_gc = slope_gc
        self.fimbh = fimbh
        self.fgcimbh = fgcimbh
        self.tdelay = tdelay
        self.scale = scale
        self.z_gc = z_gc
        self.sigma_gc = sigma_gc


    def get_galaxy_mass(self):

        ymin = (self.mgal_max / self.mstar) ** (-self.astar) * np.exp(-self.mgal_max / self.mstar)
        ymax = (self.mgal_min / self.mstar) ** (-self.astar) * np.exp(-self.mgal_min / self.mstar)
        yran, yfun = 0.0, 0.0

        while yran >= yfun:
            mgal = 10 ** (random.random() * (np.log10(self.mgal_max) - np.log10(self.mgal_min)) + np.log10(self.mgal_min))
            yran = 10 ** (random.random() * (np.log10(ymax) - np.log10(ymin)) + np.log10(ymin))
            yfun = (mgal / self.mstar) ** (-self.astar) * np.exp(-mgal / self.mstar)

        return mgal


    def get_galaxy_radius(self, mgal):

        if self.galaxy_type == 'late':
            a, b, c, m0 = 0.14, 0.39, 0.1, 3.98e10
            radius = 10 ** (np.log10(c) + a * np.log10(mgal) + (b - a) * np.log10(1. + mgal / m0))
        elif self.galaxy_type == 'early':
            a, b = 0.56, 2.88e-6
            radius = 10 ** (np.log10(b) + a * np.log10(mgal))

        return radius


    def mnsc_mgal(self, mgal):
    
        if self.galaxy_type == 'late':
            al, be, cc1, cc2, sigma = 1.001, 0.016, 2.78 * 10 ** 6.0, 3.94 * 10 ** 9.0, 0.127  # late
            all, alh, bel, beh = 0.934, 1.055, -0.045, 0.039
        elif self.galaxy_type == 'early':
            al, be, cc1, cc2, sigma = 1.363, 0.010, 2.24 * 10 ** 6.0, 1.75 * 10 ** 9.0, 0.157  # early
            all, alh, bel, beh = 1.292, 1.492, -0.050, 0.057

        mcl = cc1 * 10 ** (al * np.log10(mgal / cc2) + be)  # NSC mass in Msun
        mcl = 2.0 * sigma * random.random() + mcl - sigma

        return mcl


    def rnsc_mgal(self, mnsc):

        if self.galaxy_type == 'late':
            al, be, cc1, cc2, sigma = 0.321, -0.011, 3.31, 3.60 * 10 ** 6.0, 0.133  # late
            all, alh, bel, beh = 0.283, 0.368, -0.042, 0.003
        elif self.galaxy_type == 'early':
            al, be, cc1, cc2, sigma = 0.347, -0.024, 6.27, 1.95 * 10 ** 6.0, 0.131  # early
            all, alh, bel, beh = 0.323, 0.371, -0.045, -0.002

        rh = cc1 * 10 ** (al * np.log10(mnsc / cc2) + be)  # half-mass radius in pc
        rh = 2.0 * sigma * random.random() + rh - sigma

        return rh


    def msmbh_mnsc(self, mgal, mnsc):

        A, B, sigma = 1.24, -12.6, 0.79

        msmbh = mnsc * 10 ** (A * np.log10(mgal) + B)
        msmbh = 2.0 * sigma * random.random() + msmbh - sigma

        return msmbh


    def get_mgc(self, mmax, mmin):

        gammaa = 1.0 - self.slope_gc

        if self.slope_gc == 1:
            mb = mmin * np.exp(random.random() * np.log(mmax / mmin))
        else:
            mb = (random.random() * (mmax ** gammaa - mmin ** gammaa) + mmin ** gammaa) ** (1.0 / gammaa)

        return mb


    def fracout_mnsc(self, mnsc):

        alpha, beta, sigma = 7.28, 0.34, 0.12
        x = np.log10(mnsc)

        fin = beta * np.tanh(x - alpha) + (1 - beta)
        finr = 2.0 * sigma * random.random() + fin - sigma
        fout= 1.0 - finr

        return fout


    def star_formation_gc(self):

        ymax, ymin = 1.0, 0.0
        yran, yfun = 0.0, 0.0
        xmin, xmax = 0.0, 6.0

        while yran >= yfun:

            xx = random.random() * (xmax - xmin) + xmin
            yran = random.random() * (ymax - ymin) + ymin
            yfun = np.exp(-(xx - self.z_gc) ** 2 / (2.0 * self.sigma_gc ** 2))

        return xx


    def rhsmbh_msigma(self, mbh):

        rh = 1.72 * np.sqrt(mbh / 4e6)

        return rh


    def smbhimbh_mergers(self):

        # ---------- initializing quantities

        amgal = []
        amnsc = []
        arnsc = []
        amsmbh = []
        arhsmbh = []
        afout = []
        amgc = []
        azform = []
        atform = []
        amimbh = []
        attot = []
        aii = []
        akkk = []

        # ---------- calculating mergers

        mnsc = 0.0
        msmbh = 0.0

        k = 0
        while mnsc < self.mnsc_min or mnsc > self.mnsc_max or msmbh < 1e5:
            k += 1
            mgal = self.get_galaxy_mass() # msun
            mnsc = self.mnsc_mgal(mgal) # msun
            msmbh = self.msmbh_mnsc(mgal, mnsc) # msun
        rnsc = self.rnsc_mgal(mnsc) # pc

        rhsmbh = self.rhsmbh_msigma(msmbh) # pc
        fout = self.fracout_mnsc(mnsc)

        tot_gcmass = fout * mnsc
        part_gcmass = 0.0
        ii = 0
        while part_gcmass < tot_gcmass:

            ii = ii + 1

            mgcmax = min(self.mgc_max, tot_gcmass) # avoid sampling a GC that is larger than the total ex-situ mass of NSC
            mgcmin = self.mgc_min
            if self.scale == 1:
                rgal = self.get_galaxy_radius(mgal) # in kpc
                mgcmin = 1e6 * (rgal/ 4.) ** 2. # scale as Fahrion+ 2021 from Milky Way values
            mgc = self.get_mgc(mgcmax, mgcmin) # msun
            part_gcmass = part_gcmass + mgc

            if random.random() < self.fgcimbh:

                zform = self.star_formation_gc()
                tform = Planck15.age(zform).value # Gyr
                mimbh = self.fimbh * mgc

                t1, t2, tmed = 0.5, 10.0, 1.0 # minimum, maximum, avaerage time
                if self.tdelay == 1: # exponential distribution
                    ttt = - 2.0 * np.log(random.random() * (np.exp(-t2 / tmed) - np.exp(-t1 / tmed))  + np.exp(-t1 / tmed))
                elif self.tdelay == 2: # 1/t distribution
                    ttt = t1 * np.exp(random.random() * np.log(t2 / t1))
                elif self.tdelay == 3: # uniform in t
                    ttt = random.random() * (t2 - t1) + t1
                ttot = ttt 

            else:

                zform = 0.0
                tform = 0.0
                mimbh = 0.0
                ttot = 0.0

            amgal.append(mgal)
            amnsc.append(mnsc)
            arnsc.append(rnsc)
            amsmbh.append(msmbh)
            arhsmbh.append(rhsmbh)
            afout.append(fout)
            amgc.append(mgc)
            azform.append(zform)
            atform.append(tform)
            amimbh.append(mimbh)
            attot.append(ttot)
            aii.append(ii)
            akkk.append(k)

        return(
            amgal,
            amnsc,
            arnsc,
            amsmbh,
            arhsmbh,
            afout,
            amgc,
            azform,
            atform,
            amimbh,
            attot,
            aii,
            akkk,
        )


    def eval(self, nsample):

        fdata_name = self.directory+'/datall.txt'
        fdata = open(fdata_name, "w")

        for k in range(nsample):
            amgal, amnsc, arnsc, amsmbh, arhsmbh, afout, amgc, azform, atform, amimbh, attot, aii, akkk = self.smbhimbh_mergers()

            for j in range(len(aii)):
                print('\t %3i \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %.3E \t %3i \t %3i' % (k, amgal[j], amnsc[j], arnsc[j], amsmbh[j], arhsmbh[j], afout[j], amgc[j], azform[j], atform[j], amimbh[j], attot[j], aii[j], akkk[j]), file=fdata)

        fdata.close()

        return None


    def __call__(self, nsample):
        ''' Compute merger rate as a function of redshift '''

        return self.eval(nsample)



