# This file is part of the SM_ET_coupling module of the MDTF code package (see LICENSE.txt)

#============================================================
# Coupling between soil moisture (SM) and evapotanspiration (ET) in summer
# Sample code to call R from python
# Code written by Alexis Berg
#
# This module calculates the correlations between SM and ET, as in Berg and Sheffield (2018), Fig.1a.
#
# Reference:
# Berg and Sheffield (2018), Soil moisture-evapotranspiration coupling in CMIP5 models: relationship with simulated climate and projections, Journal of Climate, 31(12), 4865-4878. 
#============================================================

import os
import subprocess
import time

#============================================================
# generate_ncl_plots - call a nclPlotFile via subprocess call
#============================================================
def generate_R_plots(RPlotFile):
    """generate_plots_call - call a RPlotFile via subprocess call
   
    Arguments:
    RPlotFile (string) - full path to R plotting file name
    """
    # check if the RPlotFile exists - 
    # don't exit if it does not exists just print a warning.
    try:
        pipe = subprocess.Popen([
            'Rscript --verbose --vanilla {}'.format(RPlotFile)
        ] , shell=True, stdout=subprocess.PIPE)
        output = pipe.communicate()[0].decode()
        print('R routine {0} \n {1}'.format(RPlotFile,output))            
        while pipe.poll() is None:
            time.sleep(0.5)
    except OSError as e:
        print('WARNING',e.errno,e.strerror)

    return 0

if os.path.isfile(os.environ["MRSOS_FILE"]):
    print("monthly soil moisture file found")

    print("computing soil moisture memory")


#============================================================
# Call R code here
#============================================================
    print("--------- Starting soil moisture memory generate figures (using R)----------------------------")
    if ( True ):
        generate_R_plots(os.environ["POD_HOME"]+"/soil_moisture_memory.R")
    else:
        print("WARNING: For testing purposes, skipping soil moisture memory figure generation")

    print("--------- Finished soil moisture memory generate figures----------------------------")

else:
    print("Soil moisture file NOT found, skip soil moisture memory")            
