# import libraries
import os
import sys
from statsmodels.tsa.arima.model import ARIMA
from scipy import stats, special
import numpy as np
import pandas as pd

# from datetime import datetime

args = sys.argv
# args = ["", "mac_loc", "climate_scenarios", "historic", "full", "sq", "1"]

# set working directory
# os.chdir(wd)

expName = args[1]
ver = args[2]
season = args[3]
skill = args[4]
nseeds = int(args[5])
# startcent = int(args[7])

# if season == "full":
#     exp = expName + "/full"
# elif season != "full":
#     exp = expName + "/seasonal"

os.makedirs("../input/" + expName + "/long_forecast/" + season, exist_ok=True)

# -----------------------------------------------------------------------------
# short forecast generator
# -----------------------------------------------------------------------------

# calculate the antecedent average annual net total supply
def calculate_annNTS(x):

    # input: average of 48 previous quarter-months of nts
    # output: forecast of annual average nts

    # transform annual average
    boxcox = np.log(x) - 8.856

    # AR(1) to transformed annual average and back transform
    annnts = np.exp(0.7382 * boxcox + 8.856)

    return annnts


# status quo method of defining indicator and confidence
def define_indicator_confidence(x):

    # input: forecast of annual average nts
    # ouput: indicator and confidence of forecast

    # upper and lower limits based on antecedent conditions
    up99limit = x + 189
    up50limit = x + 50
    low99limit = x - 189
    low50limit = x - 50

    # conditions for wet and dry indicators [ECCC code lines 157-158]
    dry = 6859
    wet = 7237

    # define indicator of wet (1), dry (-1), or neither (0) for supply
    if x > wet:
        indicator = 1
    elif x >= dry and x <= wet:
        indicator = 0
    else:
        indicator = -1

    # compute the confidence level
    if indicator == 1:
        if low99limit >= wet:
            confidence = 3
        elif low50limit >= wet:
            confidence = 2
        elif low50limit < wet:
            confidence = 1
        else:
            confidence = np.nan

    elif indicator == 0:
        if low99limit >= dry and up99limit <= wet:
            confidence = 3
        elif low50limit >= dry and up50limit <= wet:
            confidence = 2
        elif low50limit < dry or up50limit > wet:
            confidence = 1
        else:
            confidence = np.nan

    elif indicator == -1:
        if up99limit <= dry:
            confidence = 3
        elif up50limit <= dry:
            confidence = 2
        elif up50limit > dry:
            confidence = 1
        else:
            confidence = np.nan

    return [indicator, confidence]


# creates synthetic forecasts of some skill
def synthetic_forecast_generator(x, varscale, season):

    # input: time series of status quo residuals for a given simulation
    # output: synthetically generated residuals of specified skill

    trans = 2000
    resids = x

    # box cox transformation on residuals
    bc = stats.boxcox(resids + trans)
    l = bc[1]

    # center transformed residuals on mean
    bcResids = bc[0]
    bcResids = (bc[0] - bcResids.mean()) / bc[0].std()
    bcResidsMean = bcResids.mean()
    bcResids = bcResids - bcResidsMean

    # initialize ARIMA(2, 0, 4) and fit coefficients
    model = ARIMA(bcResids, order=(2, 0, 4))
    modelPar = model.fit()

    # simulate new residuals
    simResids = model.simulate(params=modelPar.params, nsimulations=len(bcResids))

    if season == "full":
        simResidsSkill = simResids * varscale

    elif season != "full":

        # skill of status quo
        simResidsSkill = simResids

        # create index of quarter-months to seasons
        seasonInd = np.concatenate(
            (
                ["winter"] * 8,
                ["spring"] * 12,
                ["summer"] * 12,
                ["fall"] * 12,
                ["winter"] * 4,
            )
        )
        seasonInd = list(islice(cycle(seasonInd), len(resids)))

        # re-scale seasonal skills by the desired amount
        soiInd = [i for i, x in enumerate(seasonInd) if x == season]
        simResidsSkill[soiInd] = simResidsSkill[soiInd] * varscale

    # back-transform residuals
    bcSim = simResidsSkill + bcResidsMean
    residsSim = special.inv_boxcox(bcSim, l) - trans
    return residsSim


# run through input hydrologic files
filelist = os.listdir("../input/" + expName + "/hydro")

for i in range(len(filelist)):

    fn = filelist[i].split(".txt")[0]
    print(fn)

    os.makedirs(
        "../input/"
        + expName
        + "/long_forecast/"
        + season
        + "/"
        + str(skill)
        + "/"
        + fn,
        exist_ok=True,
    )

    # load input data
    data = pd.read_table("../input/" + expName + "/hydro/" + filelist[i])
    data = data.loc[:, ["Sim", "QM", "ontNTS"]]
    data.columns = ["sm", "qm", "obsNTS"]

    # format skill to numeric if not status quo
    if skill != "sq":
        skill = float(skill)

    # number of simulations
    t = data.shape[0]

    for p in range(nseeds):

        # -----------------------------------------------------------------------------
        # status quo forecast
        # -----------------------------------------------------------------------------

        if skill != 0:

            # average previous annual nts --> forecast of annual average nts
            statusquo = data.copy()
            statusquo["pavgNTS"] = statusquo.obsNTS.rolling(
                48, min_periods=48, center=False, closed="left"
            ).mean()
            statusquo = statusquo.dropna().reset_index(drop=True)
            statusquo["sqNTS"] = statusquo.pavgNTS.apply(calculate_annNTS)

            # define indicator and confidence
            ind_conf = statusquo.sqNTS.apply(define_indicator_confidence)
            statusquo["indicator"] = ind_conf.apply(pd.Series).iloc[:, 0]
            statusquo["confidence"] = ind_conf.apply(pd.Series).iloc[:, 1]

        # stop and save status quo predictions
        if skill == "sq":
            longForecast = statusquo[["sm", "qm", "sqNTS", "indicator", "confidence"]]
            longForecast.columns = ["Sim", "QM", "forNTS", "indicator", "confidence"]

        # -----------------------------------------------------------------------------
        # perfect forecast
        # -----------------------------------------------------------------------------

        elif skill != "sq":

            # average next annual nts
            perfect = data.copy()
            perfect["perNTS"] = (
                perfect.obsNTS[::-1]
                .rolling(48, min_periods=48, center=False, closed="left")
                .mean()[::-1]
            )
            perfect = perfect.dropna().reset_index(drop=True)

            # define indicator and confidence
            ind_conf = perfect.perNTS.apply(define_indicator_confidence)
            perfect["indicator"] = ind_conf.apply(pd.Series).iloc[:, 0]
            perfect["confidence"] = ind_conf.apply(pd.Series).iloc[:, 1]

            # stop and save perfect predictions
            if skill == 0:
                longForecast = perfect[
                    ["sm", "qm", "perNTS", "indicator", "confidence"]
                ]
                longForecast.columns = [
                    "Sim",
                    "QM",
                    "forNTS",
                    "indicator",
                    "confidence",
                ]

            # -----------------------------------------------------------------------------
            # synthetically generated forecast
            # -----------------------------------------------------------------------------

            synthetic = statusquo.merge(perfect, on=["sm", "qm"])
            synthetic = synthetic.loc[:, ["sm", "qm", "perNTS", "sqNTS"]]
            synthetic["orgRes"] = synthetic["perNTS"] - synthetic["sqNTS"]
            synthetic = synthetic.dropna().reset_index(drop=True)
            synthetic["simRes"] = synthetic_forecast_generator(
                synthetic.orgRes, skill, season
            )
            synthetic["synNTS"] = synthetic["perNTS"] - synthetic["simRes"]
            synthetic = synthetic.loc[:, ["sm", "qm", "synNTS"]]

            # define indicator and confidence
            ind_conf = synthetic.synNTS.apply(define_indicator_confidence)
            synthetic["indicator"] = ind_conf.apply(pd.Series).iloc[:, 0]
            synthetic["confidence"] = ind_conf.apply(pd.Series).iloc[:, 1]

            # save synthetic forecasts
            longForecast = synthetic[["sm", "qm", "synNTS", "indicator", "confidence"]]
            longForecast.columns = ["Sim", "QM", "forNTS", "indicator", "confidence"]

        # save final output file
        longForecast.to_csv(
            "../input/"
            + expName
            + "/long_forecast/"
            + season
            + "/"
            + str(skill)
            + "/"
            + fn
            + "/skill_"
            + str(skill)
            + "_S"
            + str(p + 1)
            + ".txt",
            index=False,
            sep="\t",
        )
