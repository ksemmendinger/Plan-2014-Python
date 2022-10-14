# Plan 2014 (Python)
**Experimental version of the Plan 2014 regulation model and Lake Ontario - St. Lawrence River (LOSLR) impact models written Python**

This repo contains code for the regulation model (water levels and flows) and impact models (coastal impacts, meadow marsh area, etc.) for the LOSLR basin. The forecast modules have been removed from the regulation model script and are now standalone scripts with functionality to change the skill of the forecasts fed into Plan 2014 (see below for more detail). This code is intended to be a computationally inexpensive simulation model for research purposes only and does not contain a module that captures Board deviations during periods of extreme high or low water levels. The repo contains historic supply input from 1900 - 2021 for out-of-the-box simulations.

## Running Plan 2014:

1. **Create conda environment with dependencies in environment.yml**

    ```
    conda env create -f environment.yml
    ```
 
2. **Create input files**

    **Hydro Files:** Input hydrologic files are provided for the historic supply data from 1900 - 2021 (`input/historic/hydro`). In order to run Plan 2014, you will need time series of:
    
    | Variable Name | Description |
    | --- | --- |
    | Sim | Simulation time step |
    | Year | Simulation year |
    | Month | Simulation month |
    | QM | Simulation quarter-month |
    | ontNTS | True Ontario net total supply |
    | ontNBS | True Ontario net basin supply |
    | erieOut | True Lake Erie outflows |
    | stlouisontOut | True Lac St. Louis - Lake Ontario flows [abstraction of Ottawa River flows] |
    | iceInd | Ice indicator [0 = no ice, 1 = formed/stable ice, 2 = unstable/forming ice] |
    | lsdamR | Roughness coefficient at Long Sault Dam |
    | saundershwR | Roughness coefficient at the headwaters of Moses-Saunders Dam |
    | ptclaireR | Roughness coefficient at Pointe-Claire |
    | ogdensburgR | Roughness coefficient at Ogdensburg |
    | cardinalR | Roughness coefficient at Cardinal |
    | iroquoishwR | Roughness coefficient at the headwaters of Iroquois Dam |
    | iroquoistwR | Roughness coefficient at the tailwaters of Iroquois Dam |
    | morrisburgR | Roughness coefficient at Morrisburg |
    | saunderstwR | Roughness coefficient at the tailwaters of Moses-Saunders Dam |
    | cornwallR | Roughness coefficient at Cornwall |
    | summerstownR | Roughness coefficient at Summerstown |
    | jetty1R | Roughness coefficient at Jetty 1 |
    | varennesR | Roughness coefficient at Varennes |
    | sorelR | Roughness coefficient at Sorel |
    | stpierreR	 | 	Roughness coefficient at Saint-Pierre |
    | threeriversR | Roughness coefficient at Trois-Rivières |
    | batiscanR | Roughness coefficient at Batiscan
    
    <br> 
    
    **Spinup Files:** You will also need a year (48 QMs) of "spin up" data (`input/historic/spin_up`). This allows the long-term average water level to be calculcated as well as the previous QM32 flow. This file should include:
    
    | Variable Name | Description |
    | --- | --- |
    | ontLevel | Ending Plan 2014 simulated water level on Lake Ontario |
    | ontFlow | Ending Plan 2014 simulated release from Moses-Saunders dam |
    
    <br>
    
3. **Generate time series of supply index (long-range forecast)**

    This code includes functionality for varying the forecast skill seasonally of the long-term supply index. When you run the script, you need to specify the:
    
    1. Input directory name
    2. Forecast skill: You may specify "sq" for the status quo AR(1) forecast model or specify a new forecast skill with a scalar value [0, 1] (where 0 = a perfect forecast and 1 = a synthetically generated forecast of the status quo skill).
    3. Season of skill: You may specify a certain time of year to vary the skill specified in #2. To change the skill over the time series, specify "full", else choose a season ("spring", "summer", "fall", "winter").
    4. Number of seeds: For a "sq" and "0" skill, enter 1. Else, specify the number of synthetic forecasts you wish to generate.
    
    <br>
    
    ```
    python longForecast.py historic full sq 1
    ```
    
    <br>
    
    The output file will have a time series with the following columns:
    
    | Variable Name | Description |
    | --- | --- |
    | Sim | Simulation time step |
    | Year | Simulation year |
    | Month | Simulation month |
    | QM | Simulation quarter-month |
    | forNTS | Forecast annual average Ontario net total supply over the next 48 quarter-months from long-term forecast |
    | indicator | Whether forNTS is wet (1), dry (-1), or neither (0) |
    | confidence | Confidence in how wet or dry forNTS is [1 = not confident, 2= average confidence, 3 = very confident] |
    
    <br>

4. **Generate time series of 4 quarter-monthly supply forecasts (short-range forecast)**

    This code includes functionality for changing the forecast skill from the status quo to a perfect forecast for Lake Ontario NBS, Erie outflows, and SLON flow short-range forecasts. When you run the script, you need to specify the:
    
    1. Input directory name
    2. Forecast skill: You may specify "sq" for the status quo AR(p) forecast models or 0 (a perfect forecast).
    
    <br>
    
    ```
    python shortForecast.py historic sq
    ```
    
    <br>
    
    The output file will have a time series with the following columns:

    | Variable Name | Description |
    | --- | --- |
    | Sim | Simulation time step |
    | Year | Simulation year |
    | Month | Simulation month |
    | QM | Simulation quarter-month |
    | ontNBS_QM1 | First (of four) quarter-month forecast of Ontario net basin supply from short-term forecast |
    | ontNBS_QM2 | Second (of four) quarter-month forecast of Ontario net basin supply from short-term forecast |
    | ontNBS_QM3 | Third (of four) quarter-month forecast of Ontario net basin supply from short-term forecast |
    | ontNBS_QM4 | Fourth (of four) quarter-month forecast of Ontario net basin supply from short-term forecast |
    | erieOut_QM1 | First (of four) quarter-month forecast of Lake Erie outflows from short-term forecast |
    | erieOut_QM2 | Second (of four) quarter-month forecast of Lake Erie outflows from short-term forecast |	
    | erieOut_QM3 | Third (of four) quarter-month forecast of Lake Erie outflows from short-term forecast |
    | erieOut_QM4 | Fourth (of four) quarter-month forecast of Lake Erie outflows from short-term forecast |
    | ontNTS_QM1 | First (of four) quarter-month forecast of Ontario net total supply (Ontario net basin supply + Lake Erie outflows) from short-term forecast |
    | ontNTS_QM2 | Second (of four) quarter-month forecast of Ontario net total supply (Ontario net basin supply + Lake Erie outflows) from short-term forecast |
    | ontNTS_QM3 | Third (of four) quarter-month forecast of Ontario net total supply (Ontario net basin supply + Lake Erie outflows) from short-term forecast |
    | ontNTS_QM4 | Fourth (of four) quarter-month forecast of Ontario net total supply (Ontario net basin supply + Lake Erie outflows) from short-term forecast |
    | slonFlow_QM1 | First (of four) quarter-month forecast of Lac St. Louis - Lake Ontario flows [abstraction of Ottawa River flows] from short-term forecast |
    | slonFlow_QM2 | Second (of four) quarter-month forecast of Lac St. Louis - Lake Ontario flows [abstraction of Ottawa River flows] from short-term forecast |
    | slonFlow_QM3 | Third (of four) quarter-month forecast of Lac St. Louis - Lake Ontario flows [abstraction of Ottawa River flows] from short-term forecast |
    | slonFlow_QM4 | Fourth (of four) quarter-month forecast of Lac St. Louis - Lake Ontario flows [abstraction of Ottawa River flows] from short-term forecast |
    
    <br>

5. **Run Plan 2014**

    When you run the script, you need to specify the:
    
    1. Input directory name
    2. Version of the code to run (historic or stochastic)
    3. Season of long-forecast skill
    4. Long-forecast skill
    5. Number of seeds (from long-forecast generator)
    6. Start file (typically set as 0 unless a simulation errors out)
    
    <br>

    ```
    python plan2014.py historic historic full sq 1 0
    ```
    
    <br>
    
    *Note: If you changed the short-range forecast the status quo to a perfect skill, you will need to manually change the directory from "sq" to "0" on line 54.*

    ```
    sf = pd.read_table("../input/" + expName + "/short_forecast/0/" + filelist[f])
    ```
    
    <br>
    
    The output is saved to a created `output/` directory. The output time series has values for the following variables at each time step:
    
    | Variable Name | Description |
    | --- | --- |
    | Sim | Simulation time step |
    | Year | Simulation year |
    | Month | Simulation month |
    | QM | Simulation quarter-month |
    | ontLevel | Ending Plan 2014 simulated water level on Lake Ontario |
    | ontFlow | Ending Plan 2014 simulated release from Moses-Saunders dam |
    | flowRegime | Ending flow regime of Plan 2014 simulated release [RC = rule curve, I = ice limit, J = ramping constraints, etc...] |
    | kingstonLevel | Ending Plan 2014 simulated water level at Kingston |
    | odgensburgLevel | Ending Plan 2014 simulated water level at Ogdensburg |
    | alexbayLevel | Ending Plan 2014 simulated water level at Alexandria Bay |
    | cardinalLevel | Ending Plan 2014 simulated water level at Cardinal |
    | iroquoishwLevel | Ending Plan 2014 simulated water level at the headwaters of Iroquois Dam |
    | morrisburgLevel | Ending Plan 2014 simulated water level at Morrisburg |
    | longsaultLevel | Ending Plan 2014 simulated water level at Long Sault Dam |
    | saunderhwLevel | Ending Plan 2014 simulated water level at the headwaters of Moses-Saunders Dam |
    | saundersitwLevel | Ending Plan 2014 simulated water level at the tailwaters of Moses-Saunders Dam |
    | cornwallLevel | Ending Plan 2014 simulated water level at Cornwall |
    | summerstownLevel | Ending Plan 2014 simulated water level at Summerstown |
    | stlouisFlow | Ending Plan 2014 simulated flows at Lac Saint-Louis |
    | ptclaireLevel | Ending Plan 2014 simulated water level at Pointe-Claire (Lac Saint-Louis) |
    
    <br>
    
## Running Impact Models: 

This repo also contains impact models for coastal impacts (upstream and downstream), commercial navigation costs, hydropower production, meadow marsh acreage, and recreational boating costs. The impact models take in the direct output from the Plan 2014 simulation model. To run the objective functions, specify the following:

1. Input directory name
2. Version of the code to run (historic or stochastic)
3. Season of long-forecast skill
4. Long-forecast skill
5. Number of seeds (from long-forecast generator)

<br>

```
python objectiveSimulation.py historic historic full sq 1
```
