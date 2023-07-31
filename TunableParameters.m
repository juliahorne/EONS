%% ============================= EONS Model =============================== 
% Julia Horne, 2019

ConstantsParameters;                         % standard parameters, constants, and conversion values
ops = odeset('RelTol',1e-9,'AbsTol',1e-9);   % model run relative and absolute tolerances (choose your own adventure!)

% key tunable parameters 
inp.fNH3      = 1e-5;                        % fraction of atm NH3 - defines initial conditions (Sagan + Chyba, 1972)
inp.fCO2      = 100;                         % fraction of CO2 PAL (400 ppm) in modern era (Zhanle + Sleep, 2002) - defines initial conditions
inp.tNH3a     = 5/365;                       % yr; photochemical lifetime of NH3 (Byrne + Goldblatt 2014)
inp.toxi      = 1e8;                         % yr; oxidative weathering timescale for continents (from COPSE steady state reqs. == total cont. organic carbon/assumed oxidative weathering flux, Bergman 2004)
inp.tcarb     = 5e8;                         % yr; carbonate weathering (from est. continental CaCO3 reservoir/est weathering flux; Bergman 2004 and Wallman + Aliosi 2012)
inp.tsil      = 5e8;                         % yr; continental weathering timescale (from total est. continental silicate reservoir/LOSCAR est. CaSiO3 weathering flux; estimated 2.5x sloewer than carbonate weathering (LOSCAR 2012, Walker + Kasting 1992)) 
inp.meta      = 1e9;                         % yr; continental metamorphism timescale (from total continental OC / metamorphic OC flux; Wallman + Aliosi 2012 and Falkowski 2012) - 3e9 years should also coincide roughly with Rodinia orogenies (1-1.3 Ga)
inp.sub       = 1e8;                         % yr; lifetime of ocean crust == assumed crustal spreading timescale = subduction timescale = volcanic outgassing timescale 
inp.Ssink     = 50/365;                      % yr; POC sinking time out of surface ocean (~ 3.5 m/day; Kriest + Oschlies, 2008; Fischer + Karavas 2009 ~ 100 m/day)
inp.Dsink     = 5000/365;                    % yr; POC sinking time out of deep ocean, assuming some mass loss from degredation and, of course, much deeper 
inp.rfix      = 5;                           % fixing:assim timescale ratio (~ 5-10x longer than assimilation, less energetically favorable (Fennel et al, 2005))
inp.dSi       = v.oc.dSi/20;                 % mol/m3; tuned estimate of modern porewater dissolved [SiO2] (upper estimate <0.1 mM, Isson + Planavsky, 2018)
inp.fman      = 7;                           % multiplier for initial mantle outflux (modern ~ 3e11 mol/yr; Holland, 2002)

% define the times when model resets between pre/post photosynthesis
% evolution- model start time 0 yr == 4 Ga
% Model has to reset when photosynthesis starts up or it crashes
inp.S1        = v.td.initphoto;              % 400 million years into the run == the start of Archean eon
inp.S2        = 4e9 - inp.S1;                % duration of the rest of the run is defined by reset time