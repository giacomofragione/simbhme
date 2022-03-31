# simbhme

#### Mergers of supermassive black holes and intermediate-mass black holes


### Citation

If you use this software in a scientific publication, I kindly ask you to cite: https://ui.adsabs.harvard.edu/abs/2022arXiv220205618F/abstract

This code is developed and maintained by [Giacomo Fragione](https://giacomofragione90.wixsite.com/giacomofragione).
To report bugs, please open an issue on GitHub. If you want to contact me, send an email to `giacomo.fragione90@gmail.com`.

### Results

`simbhme` has been used in the following papers:
- Fragione (2022) [https://ui.adsabs.harvard.edu/abs/2022arXiv220205618F/abstract](https://ui.adsabs.harvard.edu/abs/2022arXiv220205618F/abstract)

### Usage

The default usage is

    import simbhme as sim
    p = sim.smbhimbh_mergers()
    p(nsample)

#### Parameters

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

#### Returns

A file `datall.txt` wiith columns

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

