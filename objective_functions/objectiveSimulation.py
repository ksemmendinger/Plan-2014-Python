# import libraries
import os
import sys
import numpy as np
import pandas as pd
from glob import glob
from datetime import datetime

args = sys.argv
# args = ["", "mac_loc", "climate_scenarios", "historic", "full", "sq", "1"]

# set working directory
# if args[1] == "mac_loc":
#     wd = "/Users/kylasemmendinger/Box/Plan_2014/objective_functions/assessment"
# os.chdir(wd)

expName = args[1]
v = args[2]
season = args[3]
skill = args[4]
nseeds = int(args[5])

# import objective functions
import objectiveFunctions

# list output files for simulation
PATH = "../simulation_model/output/" + expName + "/" + season + "/" + skill
filelist = [y for x in os.walk(PATH) for y in glob(os.path.join(x[0], "*.txt"))]

for fn in filelist:

    print(fn)

    # load plan 2014 output
    data = pd.read_csv(fn, sep="\t")

    # format output
    dataTS = data.dropna().reset_index(drop=True)
    dataTS = {x: dataTS[x].values for x in dataTS}

    (
        coastal,
        commNav,
        hydro,
        mMarsh,
        recBoat,
    ) = objectiveFunctions.objectiveSimulation(dataTS)

    output = (
        data.merge(coastal, on=["Sim", "Year", "Month", "QM"], how="left")
        .merge(commNav, on=["Sim", "Year", "Month", "QM"], how="left")
        .merge(hydro, on=["Sim", "Year", "Month", "QM"], how="left")
        .merge(mMarsh, on=["Year", "QM"], how="left")
        .merge(recBoat, on=["Sim", "Year", "Month", "QM"], how="left")
    )

    dirName = (
        "output/" + expName + "/" + season + "/" + str(skill) + "/" + fn.split("/")[-2]
    )
    # make directory for output
    os.makedirs(
        dirName,
        exist_ok=True,
    )

    output.to_csv(dirName + "/" + fn.split("/")[-1], index=False, sep="\t")
