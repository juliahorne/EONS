# EONS
## Earth Oxygenation and Natural Systematics: A biogeochemical model of Earth's surface evolution spanning 4.5 billion years

*Version 1.0 - 2023*
 
### By Julia Horne and Colin Goldblatt

*School of Earth and Ocean Sciences, University of Victoria*
#

 ## GENERAL INFO
 The EONS model uses forward modelling techniques to produced a theoretical history of Earth's surface evolution. This entails using initial conditions that roughly describe the Earth's surface at the end of the Hadean eon and movement of species (CO<sub>2</sub>, NH<sub>3</sub>, biomass, etc.) between reservoirs (atmosphere, ocean, continents, mantle, etc.) through time to arrive at a preindustrial contemporary state. 
 
The model uses a system of ordinary differential equations (ODEs) describing movement of and reactions between several species in the Earth's surface reservoirs. These ODEs are evaluated through time by MATLAB's ode15s solver for stiff systems. This model was developed in MATLAB versions 2019A and 2021A. Users running with earlier or much later versions may experience some issues with in-built function calls. 

 ## COMPONENTS
 The entire model is run in the ```Run_EONS.m``` master script; this produces a nominal run as described in the published (TBD) paper, as well as the mantle sensitivity testing and oxygen history comparison plots therein. Functions for reproducing the plots in that paper are included at the end of each section in the master script; upon downloading the .zip file and adding the necessary folders to their path, the user should be able to run only that script to return all figures from the paper.
 
 Everthing included in the .zip file falls into four categories:
 
<details>
<summary> Functions comprising the model</summary>
   Theses functions are those called by the ODE solver. 
 
   | Name | Purpose | Output structure |
   |-----:|-----------| -----------|
   |```Flux_BGC.m ```        | biogeochemical fluxes                                            | ```flux``` |
   |```Flux_AirSea.m```      | air-sea gase exchange fluxes                                     | ```gasex```|
   |```Flux_Spec.m```        | chemical speciation (carbonate and reduced nitrogen equilibria)  | ```conc``` |
   |```Flux_Mix.m```         | ocean mixing (overturn) and dissolved species diffusion          | ```mix```  |
   |```Flux_Temp.m```        | temperature and radiative forcing                                | ```rf```, ```Fs```, ```Tq``` |
   |```Forcings.m ```        | time dependent imposed model forcings                            | ```tdep``` |
   |```Timescales.m ```      | timescales for all model fluxes                                  | ```tau``` |
   |```ODEs.m ```            | system of ordinary differential equations                        | ```dy```, ```dt``` |
   |```InitialConditions.m```| initial reservoir sizes                                          | ```indx```, ```y0``` |
   
 </details>

<details>
<summary> Scripts containing key information/model inputs</summary>
 These are scripts that hold constants and parameters used in the model; they output structures that are then transfered between the model functions such that we avoid using global variables. Literature reference outputs include reservoir (rx = ranges, rxr = preferred values) and flux (rflux = ranges, fluxx = preferred values) estimates for the modern Earth system; these are called by plotting functions.
   
   | Name | Purpose | Output structure |
   |-----:|-----------| -----------|
   |```TunableParameters.m ```  |all tunable model parameters                              | ```inp``` |
   |```ConstantParameters.m ``` | all constants and conversion values                      | ```v``` |
   |```LiteratureReference.m ```|compilation of literature estaimtes for fluxes/reservoirs | ```rx```,```rxr```,```rflux```,```fluxx``` |
   
 </details>

 <details>
<summary> Ancillary functions/scripts/folders called by the main model functions </summary>
  These are used by main model functions to evaluate different relationships (ie. weathering sensitivities to CO2 and temperature) or to use for parameterizing climate effects from greenhouse gas partial pressures (ie. TempParam.m and RFInterp.m calculate climate constants and radiative forcings using the ByrneSI folder model output). 
  
   | Name | Purpose | Output structure |
   |-----:|-----------| -----------|
   | ```WeatheringSensitivities.m ``` | CO2 carbonate/silicate weathering sensitivities                   |``` ws ```|
   | ```VolumetricConcentrations.m``` | calculate dissolved species concentrations in mol/m<sup>3</sup>   | ```c``` |
   | ```Limitations.m ```             | half-saturation limitations for biogeochemical reactions          |``` lim``` |
   | ```RFInterp.m ```                | interpolate greenhouse gases onto radiative transfer model output | ```rf``` |
   | ```TempParam.m ```               | temperature constant parameterization                             | ```a1```, ```b1```, ```q ```|
   | ```UnpackOutput.m ```            | reruns model fluxes with dt, dy ode outputs and calculates fluxes, reservoirs, etc.| ```t```,```r```,```flux```, etc. |
   | ```TotalOceanFluxes.m```         | sum specified fluxes across all ocean boxes                       |``` tf``` |
   | ```TotalOceanReservoirs.m```     | sum specified reservoirs across all ocean boxes                   | ```tr ```|
   </details>

 <details>
  
 <summary>Plotting functions</summary>
   These are the functions and scripts producing the plots shown in the EONS paper results section. 
  
   | Name | Purpose |
   |-----:|-----------|
   | ```Plot_NominalRun.m```          | plot all results figures showing reservoirs, fluxes through time |
   | ```Plot_MantleTest.m ```         | plot tests for different treatments of mantle reductant outgassing forcing |
   | ```Plot_OxygenHistory.m ```      | plot modelled oxygen curve against compiled literature oxygen proxies |
   | ```Plot_NutrientLimitations.m``` | plot changing N and P limitations on biosphere and C:P ratio of organic matter |
   | ```PrintPDFToFolder.m ```        | generate a PDF out of the current figure and save it to a folder on the current path |  |
 </details>

  <details>
 <summary>Evaluation functions</summary>
   These are extra functions that allow the user to check that the model is functioning properly, in particular that the model is conserving mass. 
   
   | Name | Purpose | Output structure |
   |-----:|-----------| -----------|
   | ```SumAllSpecies.m ```   | total all species, elements through time in a model run (mass tracking) | ```totalres```, ```specres``` |
   | ```MassConservation.m ```| using SumAllSpecies.m, plot species/element reservoirs and net change   |  |
   | ```DetectElement.m```    | called by SumAllSpecies.m, see if element X exists in species Y         |  |
   
  </details>

## ORGANIZATION
 The model uses structures to organize outputs and inputs. 
   * Reservoirs ```r``` are organized by: ```r.(box).(species)``` - ie. atmospheric carbon dioxide reservoir is ```r.a.CO2```
   * Fluxes are organized: ```flux.(name).(species).(box)``` unless the flux is only within one reservoir, then it's ```flux.(name).(species)```; if a flux is only for one or two species, then it's ```flux.(name)``` only.
 
 Standards, constants, and conversion parameters are kept within a single structure ```v```, with sub-categories:
 
   * ```conv```  - conversion between units of measure
   * ```const``` - general constants (ie. C:N:P Redfield ratios, molar masses, Boltzman constant, etc)
   * ```ea```    - earth system constants (ie. crustal mass, radius, etc)
   * ```atm```   - atmosphere constants (modern surface pressure, thickenss, diffusion constants, etc)
   * ```oc```    - ocean constants (depth, salinity, mass, etc)
   * ```sed```   - sediment constants (reactive layer depth, diffusivity for species, etc)
   * ```f```     - constant fractions (subduction, area of shelf sediments, etc)
   * ```spx```   - speciation constants (for equilibrium speciation)
   * ```td```    - time dependent forcings constants (initialization times, duration of transitions, etc)
   * ```Ki```,```Ks``` - inhibition/sensitivity thresholds for limitations

Tunable parameters are included in the ```inp``` structure, and can be modified by the user in sensitivity tests like the example mantle reductant influx sensitivity test detailed in the paper and included in section 2 of the master script. 

Concentrations found from evaluating carbonate or nitrogen speciation states in ```Flux_Spec.m``` are returned in the ```conc``` structure, in units of mol/kg; this function computes the concentrations of implicit species, including those not explicitly tracked in the model reservoirs (ie. DIC is an explicitly tracked species, but includes the implicit species CO<sub>2</sub>, CO<sub>3</sub><sup>2-</sup>, and HCO<sub>3</sub><sup>-</sup>). Concentrations of any dissolved species in the ocean, implicit or explicit, are calculated in the ```VolumetricConcentrations.m``` function in units of mol/m<sup>3</sup> and returned in the ```c``` structure. Both ```conc``` and ```c``` use the same organization as the reservoir ```r``` structure.
     
Functions beginning with ```plot_``` denote plotting assistants; these help produce output plots consistent with the EONS publication. 

## STEP BY STEP INSTRUCTIONS FOR STARTUP
1. Download the .zip file for the EONS model repository
2. Open/download MATLAB version 2019a - 2023b (other versions may need reviewing for outdated function usage)
3. Running the ```ConstantsParameters.m``` script should add all the files and folders to your path (with the exception of ```arrow.m``` function, which is added when ```Plot_OxygenHistory.m``` is called); this script is called by the master script. Model functions/scripts/folders should appear in the *Current Folder* window; ensure that all included subfolders are highlighted. To add the EONS folder and all subfolders to your path:
   - option 1: type
     ```
     addpath ~/path/to/EONS/
     ```
     into the command window. You will additionally have to add ```ByrneSI``` folder to the path with a second command
     ```
     addpath ~/path/to/EONS/ByrneSI
     ```
   - option 2: navigate via MATLAB tabs; select
     
     *Home* >> *Environment* >> *Set Path* icon (two folders)
    
     then select the local folder that hosts the EONS model and hit *Save*

 If for some reason not all folders are not highlighted (again, excepting ```arrow.m``` for the moment), right click on the problematic one and select in the dropdown menu:

  *Add to path* >> *Selected folders and subfolders* 

4. Open ```Run_EONS.m``` and run the desired sections
   - highlight individual sections (separated by ```%%``` commented lines) by clicking on the box and use command+enter to evaluate, or in the *Editor* tab select *Run Section*
   - run the entire script at once by using *Run* in the *Editor* tab

Note: a single model run takes **15-30 minutes on average** to evaluate, and the mantle sensitivity tests include 20 such runs evaluated using MATLAB's parallel processing functionality. We recommend running each section individually.

Each section of the master script includes plotting function calls and saves the run output as a file. Saved .mat files will end up in the path folder unless otherwise designated by the user.  Section 1 produces a nominal run from Hadean to modern conditions and plots showing system evolution, oxygen history versus proxy record, and nutrient limitations. Section 2 runs the mantle reductant sensitivity test and plots the results in a single figure. 

See the attached flow chart for a visualization of how functions and scripts relate within the EONS model.

![EONS_FlowChart.pdf](https://github.com/juliahorne/EONS/files/12053244/EONS_FlowChart.pdf)
