# missing_points_phase
Replication code for ISMRM 2021 #1184

Robust and Computationally Efficient Missing Point and Phase Estimation for Zero Echo Time (ZTE) Sequences

Curtis Corum*1,2, Abdul Haseeb Ahmed2, Mathews Jacob2, Vincent Magnotta2 and Stanley Kruger2

  1Champaign Imaging LLC, Shoreview, MN, USA

  2University of Iowa, Iowa City, IA, USA

  *Curtis Corum of Champaign Imaging LLC is developing products related to this research

<a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-nd/4.0/">Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License</a>.


So far tested only in linux environment, specifically Ubuntu 18.04

    Requires: matlab

Because there are submodules please use '--recursive' to clone:
    
    git clone --recursive https://github.com/curtcorum/missing_points_phase
    
    git submodule init
    
    git submodule update

To build mexa64 binaries

    cd matlab/DCF_Estimation
    
    make

To do a gridding reconstruction of Silent Scan p-file data, etc:

    cd ngfn_simulation_recon

    nice matlab

Once matlab starts up:

    >testobj = ngfnRecon;
    
And select the .mat file to reconstruct via the dialog:
    
![file dialog for p-file](https://github.com/curtcorum/ngfn_simulation_recon/blob/ismrm2021/ngfnRecon_dialog.png)

Output including a log and nifti files will be created in the same directory as the source dataset.

Datasets with various amounts of noise and missing points can be downloaded from:

https://www.dropbox.com/sh/16fctr9otvzbj3s/AACqbi1CINd2NuMTi72tg8daa?dl=0



