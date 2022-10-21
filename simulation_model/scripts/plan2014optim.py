# import libraries
import os
import sys
import numpy as np
import pandas as pd
from glob import glob
from datetime import datetime

args = sys.argv
# args = ["", "historic", "historic"]

expName = args[1]
v = args[2]

# load in objective function evaluations
# sys.path.insert(1, "../../objectiveFunctions")
# import objectiveFunctions

# -----------------------------------------------------------------------------
# simulation functions
# -----------------------------------------------------------------------------


def engRound(value, k):

    if k < 0:
        n = abs(k)
    else:
        n = k - 1

    n = k
    answer = np.round(value, n)

    return answer


def rcSimulation(data, v, *vars):

    # take in new decision variables
    dv = vars
    dv_ca = dv[0]
    dv_cac = dv[1]
    dv_cb = dv[2]
    dv_pa = dv[3]
    dv_pb = dv[4]
    dv_50 = dv[5]
    dv_99 = dv[6]
    dv_t = dv[7]
    dv_max = dv[8]
    dv_min = dv[9]
    dv_w = dv[10]
    dv_d = dv[11]

    # initialize columns for output
    data["flowRegime"] = str(np.nan)
    data["stlouisFlow"] = np.nan
    data["kingstonLevel"] = np.nan
    data["alexbayLevel"] = np.nan
    data["brockvilleLevel"] = np.nan
    data["ogdensburgLevel"] = np.nan
    data["cardinalLevel"] = np.nan
    data["iroquoishwLevel"] = np.nan
    data["iroquoistwLevel"] = np.nan
    data["morrisburgLevel"] = np.nan
    data["longsaultLevel"] = np.nan
    data["saundershwLevel"] = np.nan
    data["saunderstwLevel"] = np.nan
    data["cornwallLevel"] = np.nan
    data["summerstownLevel"] = np.nan
    data["lerybeauharnoisLevel"] = np.nan
    data["ptclaireLevel"] = np.nan
    data["jetty1Level"] = np.nan
    data["stlambertLevel"] = np.nan
    data["varennesLevel"] = np.nan
    data["sorelLevel"] = np.nan
    data["lacstpierreLevel"] = np.nan
    data["maskinongeLevel"] = np.nan
    data["troisrivieresLevel"] = np.nan
    data["batiscanLevel"] = np.nan
    data["freshetIndicator"] = np.nan

    # number of rows (time steps) in data frame
    timesteps = data.shape[0]

    # convert to dictionary for faster lookup
    data = {x: data[x].values for x in data}

    # set start iteration at 49 to allow for one year of spin up
    s = 48

    # 2970 cms-quarters is the conversion factor for converting flows to levels
    conv = 2970

    # -----------------------------------------------------------------------------
    # plan 2014 simulation
    # -----------------------------------------------------------------------------

    sim_st = datetime.now()

    for t in range(s, timesteps):

        # quarter month
        qm = data["QM"][t]

        # -------------------------------------------------------------------------
        # starting values for time step, t
        # -------------------------------------------------------------------------

        # ontario water level
        ontLevelStart = data["ontLevel"][t - 1]

        # kingston water level
        kingLevelStart = ontLevelStart - 0.03

        # average level of previous 48 quarter-months
        annavgLevel = np.mean(data["ontLevel"][(t - 48) : (t)])

        # moses-saunders release
        ontFlowPrev = data["ontFlow"][t - 1]

        # ice status
        iceIndPrev = data["iceInd"][t - 1]

        # -------------------------------------------------------------------------
        # short-term supply forecasts over next 4 quarter-months
        # -------------------------------------------------------------------------

        # ontario net total supply (ontario nbs + erie outflows)
        sfSupplyNTS = [
            data["ontNTS_QM1"][t],
            data["ontNTS_QM2"][t],
            data["ontNTS_QM3"][t],
            data["ontNTS_QM4"][t],
        ]
        # lac st. louis flows - lake ontario flows (ottawa river flows)
        sfSupplySLON = [
            data["slonFlow_QM1"][t],
            data["slonFlow_QM2"][t],
            data["slonFlow_QM3"][t],
            data["slonFlow_QM4"][t],
        ]

        # -------------------------------------------------------------------------
        # long-term supply forecasts
        # -------------------------------------------------------------------------

        # ontario basin supply
        lfSupply = data["forNTS"][t]
        # lfCon = data["confidence"][t]
        # lfInd = data["indicator"][t]

        # -------------------------------------------------------------------------
        # long forecast generator
        # -------------------------------------------------------------------------

        # upper and lower limits based on antecedent conditions
        up99limit = lfSupply + dv_99
        up50limit = lfSupply + dv_50
        low99limit = lfSupply - dv_99
        low50limit = lfSupply - dv_50

        # conditions for wet and dry indicators
        dry = dv_d
        wet = dv_w

        # define indicator of wet (1), dry (-1), or neither (0) for supply
        if lfSupply > wet:
            lfInd = 1
        elif lfSupply >= dry and lfSupply <= wet:
            lfInd = 0
        else:
            lfInd = -1

        # compute the confidence level
        if lfInd == 1:
            if low99limit >= wet:
                lfCon = 3
            elif low50limit >= wet:
                lfCon = 2
            elif low50limit < wet:
                lfCon = 1
            else:
                lfCon = np.nan

        if lfInd == 0:
            if low99limit >= dry and up99limit <= wet:
                lfCon = 3
            elif low50limit >= dry and up50limit <= wet:
                lfCon = 2
            elif low50limit < dry or up50limit > wet:
                lfCon = 1
            else:
                lfCon = np.nan

        if lfInd == -1:
            if up99limit <= dry:
                lfCon = 3
            elif up50limit <= dry:
                lfCon = 2
            elif up50limit > dry:
                lfCon = 1
            else:
                lfCon = np.nan

        data["indicator"][t] = lfInd
        data["confidence"][t] = lfCon

        # -------------------------------------------------------------------------
        # state indicators
        # -------------------------------------------------------------------------

        # ice status
        iceInd = data["iceInd"][t]

        # roughness coefficients
        # r = data.filter(like="R").loc[t, :]

        # true versus forecasted slon
        foreInd = 0
        if foreInd == 0:
            slonFlow = data["slonFlow_QM1"][t]
        elif foreInd == 1:
            slonFlow = data["stlouisontOut"][t]

        if v == "stochastic":
            # needed for CO limit calculation
            stfrancoisFlow = data["lacstfrancoisFlow"][t]
            beauharnoisMaxFlow = data["beauharnoisMaxFlow"][t]

        # true nts
        obsontNTS = data["ontNTS"][t]

        year = data["Year"][t]

        # flow, level, and flag if september levels are dangerously high
        if qm <= 32:

            # take the flow and level from the previous year
            qm32Flow = data["ontFlow"][data["Year"] == year - 1][32 - 1]
            qm32Level = data["ontLevel"][data["Year"] == year - 1][32 - 1]

            flowflag = 0

        elif qm > 32:

            # take the flow and level from the current year
            qm32Flow = data["ontFlow"][data["Year"] == year][32 - 1]
            qm32Level = data["ontLevel"][data["Year"] == year][32 - 1]

            if qm32Level > 74.8:
                flowflag = 1
            else:
                flowflag = 0

        # -------------------------------------------------------------------------
        # rule curve release regime
        # -------------------------------------------------------------------------

        # calculate rule curve release for each forecasted quarter-month (1 - 4)
        nforecasts = 4
        startLev = []
        startLev.append(ontLevelStart)
        endLev = []
        sfFlow = []
        sfpreprojFlow = []
        sfRegime = []

        for k in range(nforecasts):

            # function of levels and long-term forecast of supplies
            slope = 55.5823

            # set indicators
            ice = 0
            adj = 0.0014 * (2010 - 1985)
            epsolon = 0.0001

            # while loop and break variables
            flg = 1
            ct = 0
            lastflow = 0

            while flg == 1:

                # only exits loop once a convergence threshold (epsolon) is met or 10
                # iterations exceeded. adjust the preproject relationship by how much the
                # long-term supply forecast varies from average

                # pre-project flows
                preproj = slope * (ontLevelStart - adj - 69.474) ** 1.5

                # above average supplies
                if lfSupply >= dv_t:

                    # set c1 coefficients based on how confident forecast is in wet
                    if lfInd == 1 and lfCon == 3:
                        c1 = dv_cac
                    else:
                        c1 = dv_ca

                    # rule curve release
                    flow = preproj + ((lfSupply - dv_t) / dv_max) ** dv_pa * c1

                    # set rc flow regime
                    if lfInd == 1 and lfCon == 3:
                        sy = "RC1+"
                    else:
                        sy = "RC1"

                # below average supplies
                if lfSupply < dv_t:

                    # set c2 coefficient
                    c2 = dv_cb

                    # rule curve release
                    flow = preproj - ((dv_t - lfSupply) / dv_min) ** dv_pb * c2

                    # set rc flow regime
                    sy = "RC2"

                # adjust release for any ice
                release = round(flow - ice, 0)

                if abs(release - lastflow) <= epsolon:
                    break

                # calculate resulting water level
                wl1 = ontLevelStart + (sfSupplyNTS[k] / 10 - release) / conv
                wl2 = wl1
                wl1 = (ontLevelStart + wl2) * 0.5
                wl = round(wl1, 2)

                # stability check
                lastflow = release
                ct = ct + 1

                if ct == 10:
                    break

            if v == "historic":

                # try to keep ontario level up in dry periods
                if annavgLevel <= 74.6:

                    # adjust release
                    release = release - 20

                    # set flow regime
                    sy = sy + "-"

            sfFlow.append(release)
            sfpreprojFlow.append(preproj)
            sfRegime.append(sy)

            # compute water level change using forecasted supply and flow
            dif1 = (sfSupplyNTS[k] / 10 - sfFlow[k]) / conv
            dif1 = engRound(dif1, 6)
            endLev1 = startLev[k] + dif1
            endLev1 = engRound(endLev1, 2)
            endLev.append(endLev1)

            # update intial conditions
            if k < 3:
                startLev.append(endLev[k])

        # compute averaged quarter-monthly release
        ontFlow = sum(sfFlow) / nforecasts
        ontFlow = engRound(ontFlow, 0)
        # dif1 = round((sfSupplyNTS[0] / 10 - ontFlow) / conv, 6)
        # dif1 = round_half_down((sfSupplyNTS[0] / 10 - ontFlow) / conv, 6)
        # ontLevel = round_half_up(ontLevelStart + dif1, 2)
        ontRegime = sfRegime[0]

        # -----------------------------------------------------------------------------
        # limit check
        # -----------------------------------------------------------------------------

        # ---------------------------------------------------------------------------
        #
        # R+ limit - dangerously high levels
        #
        # from qm 32 september check comments in ECCC fortran code: if sep 1 lake
        # levels are dangerously high (above 75.0), begin adjusting rule curve flow
        # to target 74.8 by beginning of qm 47 and sustain through qm 48. reassess
        # each qm and modify the adjustment
        #
        # ---------------------------------------------------------------------------

        if qm >= 32 and flowflag == 1 and ontLevelStart > 74.80:

            if qm <= 46:
                flowadj = ((ontLevelStart - 74.80) * conv) / (46 - qm + 1)
            else:
                flowadj = ((ontLevelStart - 74.80) * conv) / (48 - qm + 1)

            # adjust rule curve flow
            ontFlow = ontFlow + flowadj

            if qm == 33:
                ontFlow = min(ontFlow, qm32Flow)

            # adjust rule curve flow
            ontFlow = engRound(ontFlow, 0)

            # # calculate resulting water level
            # dif1 = round((sfSupplyNTS[0] / 10 - ontFlow) / conv, 6)
            # ontLevel = round(ontLevel + dif1, 2)

            # adjust rule curve flow regime
            ontRegime = "R+"

        # compute averaged quarter-monthly level still using forecast
        # dif1 = round((sfSupplyNTS[0] / 10 - ontFlow) / conv, 6)
        dif1 = engRound((sfSupplyNTS[0] / 10 - ontFlow) / conv, 6)
        ontLevel = ontLevelStart + dif1

        # -----------------------------------------------------------------------------
        #
        # I limit - ice stability
        #
        # maximum i-limit flow check. ice status of 0, 1, and 2 correspond to no ice,
        # stable ice formed, and unstable ice forming
        #
        # -----------------------------------------------------------------------------

        if iceInd == 2 or iceIndPrev == 2:
            iLimFlow = 623

        elif iceInd == 1 or (qm < 13 or qm > 47):

            # calculate release to keep long sault level above 71.8 m
            con1 = (kingLevelStart - 62.4) ** 2.2381
            con2 = ((kingLevelStart - 71.80) / data["lsdamR"][t]) ** 0.387
            qx = (22.9896 * con1 * con2) * 0.1
            iLimFlow = engRound(qx, 0)

        else:
            iLimFlow = 0

        iRegime = "I"

        # -----------------------------------------------------------------------------
        #
        # L limit - navigation and channel capacity check
        #
        # maximum l-limit flow check - primarily based on level. applied during
        # the navigation season (qm 13 - qm 47) and during non-navigation season.
        # reference table b3 in compendium report
        #
        # -----------------------------------------------------------------------------

        lFlow = 0

        # navigation season
        if qm >= 13 and qm <= 47:

            lRegime = "LN"

            if ontLevel <= 74.22:
                lFlow = 595

            elif ontLevel <= 74.34:
                lFlow = 595 + 133.3 * (ontLevel - 74.22)

            elif ontLevel <= 74.54:
                lFlow = 611 + 910 * (ontLevel - 74.34)

            elif ontLevel <= 74.70:
                lFlow = 793 + 262.5 * (ontLevel - 74.54)

            elif ontLevel <= 75.13:
                lFlow = 835 + 100 * (ontLevel - 74.70)

            elif ontLevel <= 75.44:
                lFlow = 878 + 364.5 * (ontLevel - 75.13)

            elif ontLevel <= 75.70:
                lFlow = 991

            elif ontLevel <= 76:
                lFlow = 1020

            else:
                lFlow = 1070

        # non-navigation season
        else:
            lRegime = "LM"
            lFlow = 1150

        # channel capacity check
        lFlow1 = lFlow
        if ontLevel >= 69.10:
            lFlow2 = (747.2 * (ontLevel - 69.10) ** 1.47) / 10
        else:
            lFlow2 = np.nan

        if lFlow2 < lFlow1:
            lFlow = lFlow2
            lRegime = "LC"

        lLimFlow = engRound(lFlow, 0)

        # -----------------------------------------------------------------------------
        #
        # M limit -  low level balance
        #
        # minimum m-limit flow check. minimum limit flows to balance low levels of
        # lake ontario and lac st. louis primarily for seaway navigation interests
        #
        # -----------------------------------------------------------------------------

        # slonFlow = slonFlow * 0.1

        # stochastic version of M limit to prevent too low of flows
        if v == "stochastic":

            if ontLevel <= 74.15:
                mFlow = 0.0

            else:

                # m-limit by quarter-month
                qmLimFlow = np.hstack(
                    [
                        [595] * 4,
                        [586] * 4,
                        [578] * 4,
                        [532] * 8,
                        [538] * 4,
                        [547] * 12,
                        [561] * 8,
                        [595] * 4,
                    ]
                )

                mFlow = qmLimFlow[qm - 1]

                if ontLevel < 74.20:
                    mq = 770 - 2 * (slonFlow * 0.1)

                    if mq < mFlow:
                        mFlow = mq

        # historic version of M limit to prevent too low of flows
        elif v == "historic":

            # compute crustal adjustment factor, fixed for year 2010
            adj = 0.0014 * (2010 - 1985)
            slope = 55.5823
            mq = 0

            # this part borrowed from 58DD to prevent too low St. Louis levels
            if ontLevel > 74.20:
                mq = 680 - (slonFlow * 0.1)

            elif ontLevel > 74.10 and ontLevel <= 74.20:
                mq = 650 - (slonFlow * 0.1)

            elif ontLevel > 74.00 and ontLevel <= 74.10:
                mq = 620 - (slonFlow * 0.1)

            elif ontLevel > 73.60 and ontLevel <= 74.00:
                mq = 610 - (slonFlow * 0.1)

            else:
                mq1 = 577 - (slonFlow * 0.1)
                if ontLevel >= (adj + 69.474):
                    mq2 = slope * (ontLevel - adj - 69.474) ** 1.5
                else:
                    mq2 = np.nan
                mq = min(mq1, mq2)

            mFlow = mq

        mLimFlow = engRound(mFlow, 0)
        mRegime = "M"

        # -----------------------------------------------------------------------------
        # J limit - stability check
        #
        # j-limit stability flow check. adjusts large changes between flow for coming
        # week and actual flow last week. can be min or max limit.
        #
        # -----------------------------------------------------------------------------

        # difference between rc flow and last week's actual flow
        flowdif = abs(ontFlow - ontFlowPrev)

        # flow change bounds
        jdn = 70
        jup = 70

        # increase upper j-limit if high lake ontario level and no ice
        if ontLevel > 75.20 and iceInd == 0 and iceIndPrev < 2:
            jup = 142

        # if flow difference is positive, check if maxJ applies
        if ontFlow >= ontFlowPrev:

            # upper J limit applies
            if flowdif > jup:
                jlim = ontFlowPrev + jup
                jmaxup = jlim
                jmin = 0
                jFlow = jlim

                if jup == 70:
                    jRegime = "J+"
                else:
                    jRegime = "JJ"

            # no jlim is applied, flow is RC flow
            else:
                jFlow = ontFlow
                jmaxup = 9999
                jmin = 0
                jRegime = ontRegime

        # if the flow difference is negative, check if minJ applies
        else:

            # lower J limit applies
            if flowdif > jdn:
                jlim = ontFlowPrev - 70
                jmaxup = 9999
                jmin = jlim
                jFlow = jlim
                jRegime = "J-"

            # no jlim is applied, flow is RC flow
            else:
                jFlow = ontFlow
                jmaxup = 9999
                jmin = 0
                jRegime = ontRegime

        jLimFlow = engRound(jFlow, 0)

        # -----------------------------------------------------------------------------
        # limit comparison
        # -----------------------------------------------------------------------------

        # this is either the J-limit (if applied) or the RC flow and regime
        maxLimFlow = jLimFlow
        maxLimRegime = jRegime

        # get the smallest of the maximum limits (L and I)
        maxLim = -9999

        if lLimFlow != 0:
            if maxLim < 0:
                maxLim = lLimFlow
                maxRegime = lRegime

        if iLimFlow != 0:
            if maxLim < 0 or iLimFlow < maxLim:
                maxLim = iLimFlow
                maxRegime = iRegime

        # compare rc flow or j limit with maximum limits (RC or J with L and I)
        if maxLim > 0 and maxLimFlow > maxLim:
            maxLimFlow = maxLim
            maxLimRegime = maxRegime

        # get the biggest of the minimum limits (M)
        minLimFlow = mLimFlow
        minLimRegime = mRegime

        # compare the maximum and minimum limits
        if maxLimFlow > minLimFlow:
            limFlow = maxLimFlow
            limRegime = maxLimRegime

        # if the limit reaches to this point, then take the minimum limit
        else:

            # if the M limit is greater than the smaller of the I/L limit, take the M limit
            if minLimFlow > maxLim:
                if minLimRegime == mRegime:
                    limFlow = minLimFlow
                    limRegime = minLimRegime
                else:
                    if maxLim > minLimFlow:
                        limFlow = maxLim
                        limRegime = maxRegime
                    else:
                        limFlow = minLimFlow
                        limRegime = mRegime
            else:
                limFlow = minLimFlow
                limRegime = minLimRegime

        # update ontFlow and ontRegime post limit check
        ontFlow = limFlow
        ontRegime = limRegime

        # -----------------------------------------------------------------------------
        #
        # C limit - coteau structure limit
        #
        # check to see if plan release needs to be reduced to limit flows through the hq
        # coteau control structure. assumes perfect forecast of Lake St. Francis local
        # inflows and maximum capacity of Beauharnois
        #
        # -----------------------------------------------------------------------------

        if v == "stochastic":

            coteauFlow = (ontFlow * 10) + stfrancoisFlow - beauharnoisMaxFlow

            if iceInd == 1 or iceIndPrev == 2:
                cLimFlow = 2500
            else:
                cLimFlow = 4000

            if coteauFlow > cLimFlow:
                cFlow = ((ontFlow * 10) - (coteauFlow - cLimFlow)) / 10
                cFlow = engRound(cFlow, 1)

                ontRegime = "CO"
                ontFlow = cFlow
                dif1 = engRound((sfSupplyNTS[0] / 10 - ontFlow) / conv, 6)
                ontLevel = ontLevelStart + dif1

        # -----------------------------------------------------------------------------
        #
        # F limit - downstream flooding
        #
        # f-limit levels check. calculate lac st. louis flow at levels at pt. claire
        # to determine if downstream flooding needs to be mitigated
        #
        # -----------------------------------------------------------------------------

        stlouisFlow = ontFlow * 10 + slonFlow

        # calculate pointe claire level
        ptclaireLevel = engRound(
            16.57 + ((data["ptclaireR"][t] * stlouisFlow / 604) ** 0.58), 2
        )

        # determine "action level" to apply at pointe claire
        if ontLevelStart < 75.3:
            actionlev = 22.10
            c1 = 11523.848
        elif ontLevelStart < 75.37:
            actionlev = 22.20
            c1 = 11885.486
        elif ontLevelStart < 75.50:
            actionlev = 22.33
            c1 = 12362.610
        elif ontLevelStart < 75.60:
            actionlev = 22.40
            c1 = 12622.784
        else:
            actionlev = 22.48
            c1 = 12922.906

        # estimate flow required to maintain pointe claire below action level
        if ptclaireLevel > actionlev:
            flimFlow = engRound((c1 / data["ptclaireR"][t] - slonFlow) / 10, 0)

            if flimFlow < ontFlow:
                ontFlow = flimFlow
                ontRegime = "F"

        # -------------------------------------------------------------------------
        # ontario and st. lawrence level and flow calculations
        # -------------------------------------------------------------------------

        # calculate final ontario water level after limits applied, this is true level using observed nts
        dif2 = engRound(((obsontNTS / 10) - ontFlow) / conv, 6)
        ontLevel = engRound(ontLevelStart + dif2, 2)

        # save ontario output for next iteration
        data["ontLevel"][t] = ontLevel
        data["ontFlow"][t] = ontFlow
        data["flowRegime"][t] = ontRegime

        # calculate st. lawrence levels
        stlouisFlow = ontFlow + (slonFlow / 10)
        data["stlouisFlow"][t] = stlouisFlow

        # kingston
        kingstonLevel = ontLevel - 0.03
        difLev = kingstonLevel - 62.40

        # ogdensburg
        ogdensburgLevel = engRound(
            kingstonLevel
            - data["ogdensburgR"][t]
            * pow(ontFlow / (6.328 * pow(difLev, 2.0925)), (1 / 0.4103)),
            2,
        )

        # alexandria bay
        alexbayLevel = engRound(
            kingstonLevel - 0.39 * (kingstonLevel - ogdensburgLevel), 2
        )

        # brockville
        brockvilleLevel = engRound(
            kingstonLevel - 0.815 * (kingstonLevel - ogdensburgLevel), 2
        )

        # cardinal
        cardinalLevel = engRound(
            kingstonLevel
            - data["cardinalR"][t]
            * pow(ontFlow / (1.94908 * pow(difLev, 2.3981)), (1 / 0.4169)),
            2,
        )

        # iroquois headwaters
        iroquoishwLevel = engRound(
            kingstonLevel
            - data["iroquoishwR"][t]
            * pow(ontFlow / (2.36495 * pow(difLev, 2.2886)), (1 / 0.4158)),
            2,
        )

        # saunders headwaters
        saundershwLevel1 = engRound(
            kingstonLevel
            - (
                data["saundershwR"][t]
                * pow(
                    (ontFlow * 10) / (21.603 * pow(difLev, 2.2586)),
                    (1 / 0.3749),
                )
            ),
            2,
        )

        if saundershwLevel1 > 73.87:
            if iceInd == 2:
                saundershwLevel = saundershwLevel1
            else:
                saundershwLevel = 73.783
        else:
            saundershwLevel = saundershwLevel1

        # iroquois tailwaters (dependent saunders headwaters level)
        iroquoistwLevel1 = engRound(
            kingstonLevel
            - data["iroquoistwR"][t]
            * pow(ontFlow / (2.42291 * pow(difLev, 2.2721)), (1 / 0.4118)),
            2,
        )

        iroquoistwLevel2 = engRound(
            73.78 + pow((ontFlow * 10), 1.841) / pow((73.78 - 55), 5.891), 2
        )

        if saundershwLevel == 73.783:
            iroquoistwLevel = iroquoistwLevel2
        else:
            iroquoistwLevel = iroquoistwLevel1

        # morrisburg (dependent saunders headwaters level)
        morrisburgLevel1 = engRound(
            kingstonLevel
            - (
                data["morrisburgR"][t]
                * (ontFlow / (2.39537 * (difLev ** 2.245))) ** (1 / 0.3999)
            ),
            2,
        )

        morrisburgLevel2 = engRound(
            73.78 + 6.799 * pow((ontFlow * 10), 1.913) / 811896440.84868, 2
        )

        if saundershwLevel == 73.783:
            morrisburgLevel = morrisburgLevel2
        else:
            morrisburgLevel = morrisburgLevel1

        # long sault (dependent saunders headwaters level)
        longsaultLevel1 = engRound(
            kingstonLevel
            - (
                data["lsdamR"][t]
                * (ontFlow / (2.29896 * (difLev ** 2.2381))) ** (1 / 0.3870)
            ),
            2,
        )

        longsaultLevel2 = engRound(
            73.78 + 1408000 * pow((ontFlow * 10), 2.188) / 12501578154791700, 2
        )

        if saundershwLevel == 73.783:
            longsaultLevel = engRound(longsaultLevel2, 2)
        else:
            longsaultLevel = longsaultLevel1

        # saunders tailwaters
        saunderstwLevel = engRound(
            44.50 + 0.006338 * pow((data["saunderstwR"][t] * ontFlow * 10), 0.7158),
            2,
        )

        # cornwall
        cornwallLevel = engRound(
            45.00 + 0.0756 * pow((data["cornwallR"][t] * ontFlow * 10), 0.364),
            2,
        )

        # summerstown
        summerstownLevel = engRound(
            46.10 + 0.0109 * pow((data["summerstownR"][t] * ontFlow * 10), 0.451),
            2,
        )

        # pointe-claire (lac st. louis)
        ptclaireLevel = engRound(
            16.57 + pow((data["ptclaireR"][t] * stlouisFlow / 60.4), 0.58), 2
        )

        # lery beauharnois (uses pointe-claire level)
        lerybeauharnoisLevel = ptclaireLevel

        # jetty 1 (montreal harbor)
        jetty1Level = engRound((ptclaireLevel * 1.657) + (-28.782), 2)

        # st. lambert
        stlambertLevel = engRound((ptclaireLevel * 1.583) + (-27.471), 2)

        # varennes
        varennesLevel = engRound((ptclaireLevel * 1.535) + (-26.943), 2)

        # sorel
        sorelLevel = engRound((ptclaireLevel * 1.337) + (-23.616), 2)

        # lac st. pierre
        lacstpierreLevel = engRound((ptclaireLevel * 1.366) + (-24.620), 2)

        # maskinonge (uses lac st pierre level)
        maskinongeLevel = lacstpierreLevel

        # troisrivieres
        troisrivieresLevel = engRound((ptclaireLevel * 1.337) + (-24.425), 2)

        # batiscan
        batiscanLevel = engRound((ptclaireLevel * 1.105) + (-20.269), 2)

        data["kingstonLevel"][t] = kingstonLevel
        data["alexbayLevel"][t] = alexbayLevel
        data["brockvilleLevel"][t] = brockvilleLevel
        data["ogdensburgLevel"][t] = ogdensburgLevel
        data["cardinalLevel"][t] = cardinalLevel
        data["iroquoishwLevel"][t] = iroquoishwLevel
        data["iroquoistwLevel"][t] = iroquoistwLevel
        data["morrisburgLevel"][t] = morrisburgLevel
        data["longsaultLevel"][t] = longsaultLevel
        data["saundershwLevel"][t] = saundershwLevel
        data["saunderstwLevel"][t] = saunderstwLevel
        data["cornwallLevel"][t] = cornwallLevel
        data["summerstownLevel"][t] = summerstownLevel
        data["lerybeauharnoisLevel"][t] = lerybeauharnoisLevel
        data["ptclaireLevel"][t] = ptclaireLevel
        data["jetty1Level"][t] = jetty1Level
        data["stlambertLevel"][t] = stlambertLevel
        data["varennesLevel"][t] = varennesLevel
        data["sorelLevel"][t] = sorelLevel
        data["lacstpierreLevel"][t] = lacstpierreLevel
        data["maskinongeLevel"][t] = maskinongeLevel
        data["troisrivieresLevel"][t] = troisrivieresLevel
        data["batiscanLevel"][t] = batiscanLevel

        # defines a freshet as a spring flow that starts 1.5 times bigger than the last QM flow at LSL
        # and stays a freshet until flows drop to 90% or less of the previous QM
        data["freshetIndicator"][t] = 1

    # convert back to a dataframe
    data = pd.DataFrame(data)
    return data


# -----------------------------------------------------------------------------
# run policies
# -----------------------------------------------------------------------------

# get policy formulations to simulate
policies = pd.read_csv("../policy_formulations.txt", sep="\t")

decisionVars = [
    "Above Average Rule Curve Multiplier",
    "Confident Above Average Rule Curve Multiplier",
    "Below Average Rule Curve Multiplier",
    "Above Average Rule Curve Power",
    "Below Average Rule Curve Power",
    "Long Forecast 50% Confidence Interval",
    "Long Forecast 99% Confidence Interval",
    "Rule Curve Threshold",
    "Above Average Adjustment",
    "Below Average Adjustment",
    "Wet Threshold",
    "Dry Threshold",
]

for p in range(policies.shape[0]):

    print(policies.loc[p, "Experiment"])

    # extract decision variables
    vars = policies.loc[p, decisionVars].to_list()

    # get lead-time and skill
    lt = policies.loc[p, "Lead-Time"]
    ltCode = lt.replace("-", "")
    skill = policies.loc[p, "Skill"]
    if skill == "Status Quo (LM)":
        skillCode = "sqLM"
    elif skill == "Status Quo (AR)":
        skillCode = "sqAR"
    elif skill == "Perfect":
        skillCode == "0"
    else:
        skillCode == skill

    # list hydro files in the experiment
    path = "../input/" + expName + "/hydro/*.txt"
    filelist = glob(path)

    for f in range(len(filelist)):

        # filename
        fn = filelist[f].split(".txt")[0].split("/")[-1]
        print(fn)

        # create output directory
        os.makedirs(
            "../output/" + expName + "/" + ltCode + "/" + str(skillCode) + "/" + fn,
            exist_ok=True,
        )

        # load input hydrologic data
        data = pd.read_table(filelist[f])

        # load short term forecast predictions (status quo)
        sf = pd.read_table("../input/" + expName + "/short_forecast/sq/" + fn + ".txt")

        # load long term forecast predictions (status quo)
        lf = pd.read_table(
            "../input/"
            + expName
            + "/long_forecast/"
            + ltCode
            + "/"
            + str(skillCode)
            + "/"
            + fn
            + "/skill_"
            + str(skillCode)
            + "_S1.txt"
        )

        # join input data, short forecast, and long forecast data
        data = data.merge(sf, how="outer", on=["Sim", "Year", "Month", "QM"])
        data = data.merge(
            lf.drop(columns=["Year", "Month"], errors="ignore"),
            how="outer",
            on=["Sim", "QM"],
        )

        # load spin up data for first year in the century simulation
        if expName != "climate_scenarios":
            spinup = pd.read_table("../input/" + expName + "/spin_up/" + fn + ".txt")
            data.loc[:47, ["ontLevel", "ontFlow"]] = spinup

        else:
            data.loc[:46, "ontLevel"] = 74.744
            data.loc[47, "ontLevel"] = 74.55
            data.loc[47, "ontFlow"] = 623
            data.loc[47, "ontLevel"] = 74.55
            data.loc[47, "kingstonLevel"] = 75.52
            data.loc[47, "stlouisFlow"] = 900

        output = rcSimulation(data, v, *vars)

        # save output
        output.to_csv(
            "../output/"
            + expName
            + "/"
            + ltCode
            + "/"
            + str(skillCode)
            + "/"
            + fn
            + "/S1.txt",
            sep="\t",
            index=False,
        )
