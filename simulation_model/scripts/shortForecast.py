# import libraries
import os
import sys
import pandas as pd
from glob import glob

# import numpy as np
# from datetime import datetime

args = sys.argv
# args = ["", "mac_loc", "historic", "sq"]

expName = args[1]
skill = args[2]
# ver = args[3]
# season = args[4]
# nseeds = int(args[6])
# startcent = int(args[7])

os.makedirs("../input/" + expName + "/short_forecast/" + skill, exist_ok=True)

# -----------------------------------------------------------------------------
# short forecast generator
# -----------------------------------------------------------------------------

# short forecast of lake ontario net basin supply
def shortforecast_ontNBS(ontNBS):

    # convert series to list
    ontNBS = ontNBS.to_list()

    # set counter to 26
    t = 26

    for k in range(4):
        ontNBS.append(
            (
                (ontNBS[t - 1] - 1034.56) * 0.4746
                + (ontNBS[t - 2] - 1034.56) * 0.1493
                + (ontNBS[t - 3] - 1034.56) * 0.01923
                + (ontNBS[t - 4] - 1034.56) * 0.03424
                + (ontNBS[t - 5] - 1034.56) * 0.01025
                - (ontNBS[t - 6] - 1034.56) * 0.006960
                - (ontNBS[t - 7] - 1034.56) * 0.03209
                + (ontNBS[t - 8] - 1034.56) * 0.05757
                + (ontNBS[t - 9] - 1034.56) * 0.00007484
                - (ontNBS[t - 10] - 1034.56) * 0.01390
                - (ontNBS[t - 11] - 1034.56) * 0.0007310
                - (ontNBS[t - 12] - 1034.56) * 0.01017
                - (ontNBS[t - 13] - 1034.56) * 0.02105
                - (ontNBS[t - 14] - 1034.56) * 0.009959
                - (ontNBS[t - 15] - 1034.56) * 0.01256
                + (ontNBS[t - 16] - 1034.56) * 0.01011
                - (ontNBS[t - 17] - 1034.56) * 0.008653
                - (ontNBS[t - 18] - 1034.56) * 0.03989
                - (ontNBS[t - 19] - 1034.56) * 0.01719
                - (ontNBS[t - 20] - 1034.56) * 0.01100
                - (ontNBS[t - 21] - 1034.56) * 0.02326
                - (ontNBS[t - 22] - 1034.56) * 0.02640
                - (ontNBS[t - 23] - 1034.56) * 0.01341
                - (ontNBS[t - 24] - 1034.56) * 0.006209
                - (ontNBS[t - 25] - 1034.56) * 0.01959
                - (ontNBS[t - 26] - 1034.56) * 0.04441
                + 1034.56
            )
        )

        # update counter
        t = t + 1

    # return 4 quarter-monthly forecasts
    output = ontNBS[(t - 4) : (t - 0)]
    return output


# short forecast of lake erie outflows into lake ontario
def shortforecast_erieOut(erieOut, qm):

    # convert series to list
    erieOut = erieOut.to_list()

    # seasonal adjustment factors
    seasadj = [
        -180.37,
        -272.70,
        -271.09,
        -295.12,
        -369.51,
        -375.63,
        -338.46,
        -287.24,
        -289.94,
        -213.70,
        -99.46,
        -47.30,
        77.74,
        67.11,
        149.16,
        191.50,
        313.56,
        375.59,
        387.25,
        368.72,
        404.33,
        301.34,
        231.34,
        254.57,
        258.97,
        221.39,
        181.37,
        164.27,
        124.83,
        109.91,
        95.68,
        78.11,
        31.54,
        18.16,
        -41.12,
        -37.91,
        -4.408,
        -102.42,
        -82.02,
        -149.59,
        -133.79,
        -175.46,
        -64.65,
        -83.87,
        -119.89,
        -182.01,
        -103.17,
        -85.62,
    ]

    # set counter to 21
    t = 21

    for k in range(4):

        # compute qm lag index
        lag = list(range(21))

        for j in range(21):
            lag[j] = qm - (j + 2)

            if lag[j] < 0:
                lag[j] = lag[j] + 48

        erieOut.append(
            (
                (erieOut[t - 1] - 6026.1946 - seasadj[lag[0]]) * 0.56880
                + (erieOut[t - 2] - 6026.1946 - seasadj[lag[1]]) * 0.21370
                + (erieOut[t - 3] - 6026.1946 - seasadj[lag[2]]) * 0.04085
                + (erieOut[t - 4] - 6026.1946 - seasadj[lag[3]]) * 0.01850
                + (erieOut[t - 5] - 6026.1946 - seasadj[lag[4]]) * 0.02194
                + (erieOut[t - 6] - 6026.1946 - seasadj[lag[5]]) * 0.03984
                + (erieOut[t - 7] - 6026.1946 - seasadj[lag[6]]) * 0.02599
                + (erieOut[t - 8] - 6026.1946 - seasadj[lag[7]]) * 0.03943
                - (erieOut[t - 9] - 6026.1946 - seasadj[lag[8]]) * 0.02275
                + (erieOut[t - 10] - 6026.1946 - seasadj[lag[9]]) * 0.01456
                + (erieOut[t - 11] - 6026.1946 - seasadj[lag[10]]) * 0.009643
                - (erieOut[t - 12] - 6026.1946 - seasadj[lag[11]]) * 0.007157
                + (erieOut[t - 13] - 6026.1946 - seasadj[lag[12]]) * 0.040900
                + (erieOut[t - 14] - 6026.1946 - seasadj[lag[13]]) * 0.005263
                - (erieOut[t - 15] - 6026.1946 - seasadj[lag[14]]) * 0.016580
                - (erieOut[t - 16] - 6026.1946 - seasadj[lag[15]]) * 0.025850
                - (erieOut[t - 17] - 6026.1946 - seasadj[lag[16]]) * 0.025210
                + (erieOut[t - 18] - 6026.1946 - seasadj[lag[17]]) * 0.003007
                - (erieOut[t - 19] - 6026.1946 - seasadj[lag[18]]) * 0.015910
                + (erieOut[t - 20] - 6026.1946 - seasadj[lag[19]]) * 0.016660
                + (erieOut[t - 21] - 6026.1946 - seasadj[lag[20]]) * 0.034700
                + 6026.1946
                + seasadj[qm - 1]
            )
        )

        # increase counter and qm
        t = t + 1
        qm = qm + 1
        if qm == 49:
            qm = 1

    # return 4 quarter-monthly forecasts
    output = erieOut[(t - 4) : (t - 0)]
    return output


# short forecast of slon flows
def shortforecast_slonFlow(slonFlow, qm):

    # convert series to list
    slonFlow = lslont
    slonFlow = slonFlow.to_list()

    # seasonal adjustment factors
    seasadj = [
        -279.19,
        -157.14,
        -133.82,
        -167.49,
        -125.78,
        -251.62,
        -362.73,
        -410.16,
        -211.33,
        -206.77,
        34.98,
        467.61,
        1013.00,
        1175.00,
        1337.00,
        1321.00,
        1327.00,
        1209.00,
        1021.00,
        777.78,
        549.02,
        336.17,
        163.36,
        23.78,
        -102.39,
        -228.72,
        -327.94,
        -403.19,
        -461.09,
        -520.08,
        -548.27,
        -574.02,
        -599.19,
        -603.32,
        -575.31,
        -540.46,
        -484.80,
        -455.02,
        -415.53,
        -297.60,
        -228.10,
        -172.10,
        -116.91,
        -114.51,
        -159.88,
        -136.68,
        -174.25,
        -207.71,
    ]

    # set counter to 5
    t = 5

    for k in range(4):

        # compute qm lag index
        lag = list(range(5))

        for j in range(5):
            lag[j] = qm - (j + 2)

            if lag[j] < 0:
                lag[j] = lag[j] + 48

        slonFlow.append(
            (
                (slonFlow[t - 1] - 1185.4761 - seasadj[lag[0]]) * 0.86820
                - (slonFlow[t - 2] - 1185.4761 - seasadj[lag[1]]) * 0.14340
                + (slonFlow[t - 3] - 1185.4761 - seasadj[lag[2]]) * 0.04482
                - (slonFlow[t - 4] - 1185.4761 - seasadj[lag[3]]) * 0.01166
                + (slonFlow[t - 5] - 1185.4761 - seasadj[lag[4]]) * 0.05068
                + 1185.4761
                + seasadj[qm - 1]
            )
        )

        # increase counter and qm
        t = t + 1
        qm = qm + 1
        if qm == 49:
            qm = 1

    # return 4 quarter-monthly forecasts
    output = slonFlow[(t - 4) : (t - 0)]
    return output


# run through input hydrologic files
# filelist = os.listdir("../input/" + expName + "/hydro")
path = "../input/" + expName + "/hydro/*.txt"
filelist = glob(path)

for i in range(len(filelist)):

    fn = filelist[i].split(".txt")[0].split('/')[-1]
    print(fn)

    # load input data
    data = pd.read_table(filelist[i])

    # initialize lists for output
    sf_nbs = data.loc[:, ["Sim", "Year", "Month", "QM"]]
    sf_erie = data.loc[:, ["Sim", "Year", "Month", "QM"]]
    sf_slon = data.loc[:, ["Sim", "Year", "Month", "QM"]]

    if skill == "sq":

        # run through the simulation of a century and make forecast predictions
        for x in range(data.shape[0]):

            qm = data.loc[x, "QM"]

            # ontario net basin supply
            if x > 25:
                nbs = data.loc[(x - 26) : (x - 1), "ontNBS"]
                nbs = shortforecast_ontNBS(nbs)
                sf_nbs.at[x, "ontNBS_QM1"] = nbs[0]
                sf_nbs.at[x, "ontNBS_QM2"] = nbs[1]
                sf_nbs.at[x, "ontNBS_QM3"] = nbs[2]
                sf_nbs.at[x, "ontNBS_QM4"] = nbs[3]

            # lake erie outflows
            if x > 20:
                le = data.loc[(x - 22) : (x - 1), "erieOut"]
                le = shortforecast_erieOut(le, qm)
                sf_erie.at[x, "erieOut_QM1"] = le[0]
                sf_erie.at[x, "erieOut_QM2"] = le[1]
                sf_erie.at[x, "erieOut_QM3"] = le[2]
                sf_erie.at[x, "erieOut_QM4"] = le[3]

            # st. lawrence - lac st. louis flows
            if x > 4:
                lslont = data.loc[(x - 5) : (x - 1), "stlouisontOut"]
                lslont = shortforecast_slonFlow(lslont, qm)
                sf_slon.at[x, "slonFlow_QM1"] = lslont[0]
                sf_slon.at[x, "slonFlow_QM2"] = lslont[1]
                sf_slon.at[x, "slonFlow_QM3"] = lslont[2]
                sf_slon.at[x, "slonFlow_QM4"] = lslont[3]

    elif skill == "0":

        # run through the simulation of a century and make forecast predictions
        for x in range(data.shape[0] - 4):

            sf_nbs.at[x, "ontNBS_QM1"] = data.loc[x + 1, "ontNBS"]
            sf_nbs.at[x, "ontNBS_QM2"] = data.loc[x + 2, "ontNBS"]
            sf_nbs.at[x, "ontNBS_QM3"] = data.loc[x + 3, "ontNBS"]
            sf_nbs.at[x, "ontNBS_QM4"] = data.loc[x + 4, "ontNBS"]

            sf_erie.at[x, "erieOut_QM1"] = data.loc[x + 1, "erieOut"]
            sf_erie.at[x, "erieOut_QM2"] = data.loc[x + 2, "erieOut"]
            sf_erie.at[x, "erieOut_QM3"] = data.loc[x + 3, "erieOut"]
            sf_erie.at[x, "erieOut_QM4"] = data.loc[x + 4, "erieOut"]

            sf_slon.at[x, "slonFlow_QM1"] = data.loc[x + 1, "stlouisontOut"]
            sf_slon.at[x, "slonFlow_QM2"] = data.loc[x + 2, "stlouisontOut"]
            sf_slon.at[x, "slonFlow_QM3"] = data.loc[x + 3, "stlouisontOut"]
            sf_slon.at[x, "slonFlow_QM4"] = data.loc[x + 4, "stlouisontOut"]

    shortForecast = sf_nbs.merge(sf_erie, on=["Sim", "Year", "Month", "QM"]).merge(
        sf_slon, on=["Sim", "Year", "Month", "QM"]
    )
    shortForecast["ontNTS_QM1"] = (
        shortForecast["ontNBS_QM1"] + shortForecast["erieOut_QM1"]
    )
    shortForecast["ontNTS_QM2"] = (
        shortForecast["ontNBS_QM2"] + shortForecast["erieOut_QM2"]
    )
    shortForecast["ontNTS_QM3"] = (
        shortForecast["ontNBS_QM3"] + shortForecast["erieOut_QM3"]
    )
    shortForecast["ontNTS_QM4"] = (
        shortForecast["ontNBS_QM4"] + shortForecast["erieOut_QM4"]
    )

    shortForecast.to_csv(
        "../input/"
        + expName
        + "/short_forecast/"
        + skill
        + "/"
        + fn
        + ".txt",
        index=False,
        sep="\t",
    )
