# import libraries
from multiprocessing.resource_sharer import stop
import os
import sys
from glob import glob
from statsmodels.tsa.arima.model import ARIMA
from scipy import stats, special
import numpy as np
import pandas as pd

args = sys.argv
# args = ["", "historic", "sqAR", "12month", "1"]

expName = args[1]
skill = args[2]
mode = args[3]
nseeds = int(args[4])
# startcent = int(args[7])

# if season == "full":
#     exp = expName + "/full"
# elif season != "full":
#     exp = expName + "/seasonal"

# os.makedirs("../input/" + expName + "/long_forecast/" + mode, exist_ok=True)

# -----------------------------------------------------------------------------
# long forecast generator
# -----------------------------------------------------------------------------

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

if mode == "full" or mode == "spring" or mode == "summer" or mode == "fall" or mode == "winter":
    
    # calculate the antecedent average annual net total supply
    def calculate_annNTS(x):

        # input: average of 48 previous quarter-months of nts
        # output: forecast of annual average nts

        # transform annual average
        boxcox = np.log(x) - 8.856

        # AR(1) to transformed annual average and back transform
        annnts = np.exp(0.7382 * boxcox + 8.856)

        return annnts

    # creates synthetic forecasts of some skill
    def synthetic_forecast_generator(x, varscale, mode):

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

        if mode == "full":
            simResidsSkill = simResids * varscale

        elif mode != "full":

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
            soiInd = [i for i, x in enumerate(seasonInd) if x == mode]
            simResidsSkill[soiInd] = simResidsSkill[soiInd] * varscale

        # back-transform residuals
        bcSim = simResidsSkill + bcResidsMean
        residsSim = special.inv_boxcox(bcSim, l) - trans
        return residsSim

    # run through input hydrologic files
    # filelist = os.listdir("../input/" + expName + "/hydro")
    path = "../input/" + expName + "/hydro/*.txt"
    filelist = glob(path)

    for i in range(len(filelist)):

        fn = filelist[i].split(".txt")[0].split('/')[-1]
        print(fn)

        os.makedirs(
            "../input/"
            + expName
            + "/long_forecast/"
            + mode
            + "/"
            + str(skill)
            + "/"
            + fn,
            exist_ok=True,
        )

        # load input data
        data = pd.read_table(filelist[i])
        data = data.loc[:, ["Sim", "QM", "ontNTS"]]
        data.columns = ["sm", "qm", "obsNTS"]

        # format skill to numeric if not status quo
        if skill != "sq":
            skill = float(skill)

        # number of simulations
        t = data.shape[0]

        for p in range(nseeds):

            # -----------------------------------------------------------------------------
            # status quo (ar) forecast
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
                    synthetic.orgRes, skill, mode
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
                + mode
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

elif 'month' in mode:

    # manually set nseed
    p = 0

    # lead time of interest in months
    ltName = mode
    leadTimeMonth = int(ltName.split("month")[0])
    leadTime = leadTimeMonth * 4

    # run through input hydrologic files
    # filelist = os.listdir("../input/" + expName + "/hydro")
    path = "../input/" + expName + "/hydro/*.txt"
    filelist = glob(path)

    for i in range(len(filelist)):

        fn = filelist[i].split(".txt")[0].split('/')[-1]
        print(fn)

        os.makedirs(
            "../input/"
            + expName
            + "/long_forecast/"
            + mode
            + "/"
            + str(skill)
            + "/"
            + fn,
            exist_ok=True,
        )

        # load input data
        data = pd.read_table(filelist[i])
        # data = data.loc[:, ["Sim", "QM", "ontNTS"]]
        # data.columns = ["sm", "qm", "obsNTS"]

        if skill == "sqAR":

            # arma model input: calculate rolling average of previous N nts obsevations, where N corresponds to the lead time of interest
            longForecast = data.loc[:, ["Sim", "Year", "Month", "QM", "ontNTS"]]
            longForecast["x"] = longForecast.ontNTS.rolling(
                leadTime, min_periods=leadTime, center=False, closed="left"
            ).mean()

            # box-cox and centering transformation
            cSQ = 8.856  # original mean of log(x) (used to center data around 0)
            longForecast["xBC"] = np.log1p(longForecast["x"]) - cSQ

            # AR(1)
            b0SQ = 0
            b1SQ = 0.7382
            longForecast["yhat_bc"] = b0SQ + (b1SQ * longForecast["xBC"])

            # backtransform to cms
            longForecast["forNTS"] = np.exp(longForecast["yhat_bc"] + cSQ)

            ind_conf = longForecast.forNTS.apply(define_indicator_confidence)
            longForecast["indicator"] = ind_conf.apply(pd.Series).iloc[:, 0]
            longForecast["confidence"] = ind_conf.apply(pd.Series).iloc[:, 1]

            # save forecasts
            longForecast = longForecast.loc[:, ["Sim", "QM", "forNTS", "indicator", "confidence"]]

        elif skill == "0":

            # perfect predictions
            longForecast = data.loc[:, ["Sim", "Year", "Month", "QM", "ontNTS"]]
            longForecast["forNTS"] = (
                longForecast.ontNTS[::-1]
                .rolling(leadTime, min_periods=leadTime, center=False, closed="left")
                .mean()[::-1]
            )

            # define indicator and confidence
            ind_conf = longForecast.forNTS.apply(define_indicator_confidence)
            longForecast["indicator"] = ind_conf.apply(pd.Series).iloc[:, 0]
            longForecast["confidence"] = ind_conf.apply(pd.Series).iloc[:, 1]

            # save forecasts
            longForecast = longForecast.loc[:, ["Sim", "QM", "forNTS", "indicator", "confidence"]]

        elif skill == "sqLM":

            # mean and sd from historic run
            mu_x = [7083.931076388889, 7084.358506944444,
            7085.242708333333,
            7085.806076388889,
            7086.752777777778,
            7087.233333333334,
            7087.4326388888885,
            7087.6618055555555,
            7088.020138888888,
            7088.713020833334,
            7089.365625,
            7090.061458333334,
            7090.4319444444445,
            7090.433680555555,
            7090.509375,
            7090.555381944445,
            7090.991319444444,
            7091.372222222222,
            7091.829340277778,
            7092.039756944445,
            7092.558333333333,
            7092.7737847222215,
            7092.870659722222,
            7093.133680555556,
            7093.503125,
            7093.865104166666,
            7093.975347222222,
            7094.291666666667,
            7094.551909722222,
            7094.81875,
            7094.784895833333,
            7095.106944444445,
            7095.434895833333,
            7095.477604166666,
            7095.894965277778,
            7095.9482638888885,
            7096.256076388889,
            7096.533333333334,
            7096.772048611111,
            7097.242708333333,
            7097.345486111111,
            7097.6064236111115,
            7098.0491319444445,
            7098.164930555556,
            7098.1644097222215,
            7098.489756944445,
            7098.6227430555555,
            7098.860763888889]
            sigma_x = [708.8117546400873,
 710.1782833878183,
 712.4779672020203,
 712.4225097548731,
 714.7960613595477,
 716.7279621469172,
 718.4125202138579,
 719.9988045064397,
 721.8184762184961,
 726.4742370259862,
 728.7274973726128,
 729.5812051088302,
 731.7823717066013,
 730.0706547540454,
 728.2584794223172,
 725.6918474864422,
 725.4210525627713,
 724.0152798972714,
 723.4855083853869,
 724.1061284410894,
 724.9149466685182,
 724.7702234931264,
 725.2468001767845,
 724.5855205505818,
 725.8508812767759,
 728.0348929638724,
 727.9042440543712,
 728.4543049941956,
 728.8773430452776,
 729.0170949215062,
 730.522753862059,
 731.4259381832112,
 732.8540688423561,
 732.8835265155097,
 732.5143204458807,
 732.6208506054708,
 733.4114344163272,
 733.0965442469752,
 732.9232812447638,
 732.3789872653095,
 731.8598014094731,
 730.3435847801035,
 728.9627143277668,
 727.0217835634953,
 726.3651125673965,
 727.5675282692239,
 725.0616592795445,
 723.8132746880291]
            mu_y = [7084.358506944444,
 7085.242708333333,
 7085.806076388889,
 7086.752777777778,
 7087.233333333334,
 7087.4326388888885,
 7087.6618055555555,
 7088.020138888888,
 7088.713020833334,
 7089.365625,
 7090.061458333334,
 7090.4319444444445,
 7090.433680555555,
 7090.509375,
 7090.555381944445,
 7090.991319444444,
 7091.372222222222,
 7091.829340277778,
 7092.039756944445,
 7092.558333333333,
 7092.7737847222215,
 7092.870659722222,
 7093.133680555556,
 7093.503125,
 7093.865104166666,
 7093.975347222222,
 7094.291666666667,
 7094.551909722222,
 7094.81875,
 7094.784895833333,
 7095.106944444445,
 7095.434895833333,
 7095.477604166666,
 7095.894965277778,
 7095.9482638888885,
 7096.256076388889,
 7096.533333333334,
 7096.772048611111,
 7097.242708333333,
 7097.345486111111,
 7097.6064236111115,
 7098.0491319444445,
 7098.164930555556,
 7098.1644097222215,
 7098.489756944445,
 7098.6227430555555,
 7098.860763888889,
 7099.646006944445]
            sigma_y = [710.1782833878183,
 712.4779672020203,
 712.4225097548731,
 714.7960613595477,
 716.7279621469172,
 718.4125202138579,
 719.9988045064397,
 721.8184762184961,
 726.4742370259862,
 728.7274973726128,
 729.5812051088302,
 731.7823717066013,
 730.0706547540454,
 728.2584794223172,
 725.6918474864422,
 725.4210525627713,
 724.0152798972714,
 723.4855083853869,
 724.1061284410894,
 724.9149466685182,
 724.7702234931264,
 725.2468001767845,
 724.5855205505818,
 725.8508812767759,
 728.0348929638724,
 727.9042440543712,
 728.4543049941956,
 728.8773430452776,
 729.0170949215062,
 730.522753862059,
 731.4259381832112,
 732.8540688423561,
 732.8835265155097,
 732.5143204458807,
 732.6208506054708,
 733.4114344163272,
 733.0965442469752,
 732.9232812447638,
 732.3789872653095,
 731.8598014094731,
 730.3435847801035,
 728.9627143277668,
 727.0217835634953,
 726.3651125673965,
 727.5675282692239,
 725.0616592795445,
 723.8132746880291,
 723.209371671723]
            
            stand = pd.DataFrame({
                'QM': list(range(1, 49)),
                'mu_x': mu_x, 
                'sigma_x': sigma_x,
                'mu_y': mu_y,
                'sigma_y': sigma_y})
            
            # coefficients
            if leadTimeMonth == 1:
                b0 = 0.001896091232479278
                b1 = 0.6816867525759694
            elif leadTimeMonth == 3:
                b0 = -0.014933155418542947
                b1 = 0.724392669243674
            elif leadTimeMonth == 6:
                b0 = -0.0038214897473424518
                b1 = 0.6976015067700171
            elif leadTimeMonth == 12:
                b0 = -0.0020442065242214664
                b1 = 0.7231412941202915

            # make predictions
            longForecast = data.loc[:, ["Sim", "Year", "Month", "QM", "ontNTS"]]
            longForecast["x"] = longForecast.ontNTS.rolling(
                leadTime, min_periods=leadTime, center=False, closed="left"
            ).mean()

            # standardize input data by quarter-monthly mean and standard deviation
            longForecast = (
                pd.merge(longForecast, stand, on=["QM"])
                .sort_values("QM")
                .reset_index(drop=True)
            )

            # make predictions over entire simulation period
            longForecast["xStand"] = (longForecast["x"] - longForecast["mu_x"]) / longForecast[
                "sigma_x"
            ]

            longForecast["yStand"] = b0 + (longForecast["xStand"] * b1)

            # unstandardize values
            longForecast["forNTS"] = (
                longForecast["yStand"] * longForecast["sigma_y"]
            ) + longForecast["mu_y"]

            # format output
            longForecast = longForecast.loc[:, ["Sim", "Year", "Month", "QM", "forNTS"]]
            longForecast = longForecast.sort_values(by=['Sim']).reset_index(drop=True)

            # define indicator and confidence
            ind_conf = longForecast.forNTS.apply(define_indicator_confidence)
            longForecast["indicator"] = ind_conf.apply(pd.Series).iloc[:, 0]
            longForecast["confidence"] = ind_conf.apply(pd.Series).iloc[:, 1]

        # save final output file
        longForecast.to_csv(
            "../input/"
            + expName
            + "/long_forecast/"
            + mode
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
