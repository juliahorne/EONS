%% ============================= EONS Model =============================== 
% Julia Horne, 2018

% based on model inputs, define rate constants for biogeochemical fluxes in 
% the ocean/sediment reservoirs
function [tau] = Timescales(tdep,inp,v)

    tau.woxi            = inp.toxi;                 % yr; oxidative weathering timescale for continents (from COPSE steady state reqs. == total cont. organic carbon/assumed oxidative weathering flux, Bergman 2004)
    tau.wcarb           = inp.tcarb;                % yr; carbonate weathering (from est. continental CaCO3 reservoir/est weathering flux; Bergman 2004 and Wallman + Aliosi 2012)
    tau.wsil            = inp.tsil;                 % yr; continental weathering timescale (from total est. continental silicate reservoir/LOSCAR est. CaSiO3 weathering flux; estimated 2.5x faster than silicate weathering (LOSCAR 2012, Walker + Kasting 1992))
    tau.wcarbP          = tau.wcarb ./ tdep.fungi;  % yr; phosphorus mining (faster P weathering) is time-dependent (starts after colonization of land by fungi ~ 700 Ma)
    tau.wsilP           = tau.wsil ./ tdep.fungi;   % yr; (see above)
    tau.meta            = inp.meta;                 % yr; continental metamorphism timescale (from total continental OC / metamorphic OC flux; Wallman + Aliosi 2012 and Falkowski 2012)
    tau.pholys          = 10;                       % yr; theoretical NH3 lifetime via anoxic photolysis (Sagan + Chyba, 1977 estimate)
    tau.ammox           = inp.tNH3a;                % yr; NH3 atmospheric modern lifetime (oxidizing timescale, operates 1-5 days; Byrne + Goldblatt, 2014)
    tau.haze            = 2.35;                     % yr; atmospheric residence of fractal photochemical haze, assuming 1e14 g/yr production (Wolf and Toon, 2010 SI table 1)
    tau.assim           = 0.25;                     % yr; fixed N assimilation timescale (calculated by tau = LB/F_assim ~ 0.2-0.5 yr, using Gruber 2008 assim flux converted to mol C/yr)
    tau.fix             = tau.assim .* inp.rfix;    % yr; N2 fixation takes 5-10x longer than assimilation (Fennel et al, 2005)
    tau.death           = 15/365;                   % yr; lifetime of living biomass (~ 2 weeks, from IPCC Ar5, 2014)
    tau.sink.s          = inp.Ssink ./ tdep.sink;   % yr; sinking timescale of surface ocean POC (exports  over 100m to d ocean) - with org size modifier depending on geologic time cicra Ediacaran large orgs
    tau.sink.d          = inp.Dsink ./ tdep.sink;   % yr; sinking timescale of deep ocean POC (exports over 4000m to z seds) - ^^ same mod
    tau.sink.z          = 1e3;                      % yr; residence time in deep reactive sediments (pelagic sedimentation 5m/1e6 - 200m/1e6 yr; Hess + Schacht, 2011) - we assume ~1-5e3 yrs as intermediate modern
    tau.sink.n          = 1e2;                      % yr; residence time in shallow reactive sediments (^^, but faster assuming 10x higher sediment influx)
    tau.carbsink.s      = tau.sink.s / 10;          % yr; calcareous shells fall much faster bro
    tau.carbsink.d      = tau.sink.d / 10;          % yr; (as above, concerning the export through deep ocean)
    tau.bifsink.s       = v.oc.depth.s / 2e4;       % yr; iron oxide sinking timescale to deep ocean, given euphotic zone depth and avg sinking rate ~2e4 m/yr (Thompson et al. 2019)
    tau.bifsink.d       = v.oc.depth.d / 2e4;       % yr; (as above, but timescale from deep ocean to deep seds)
    tau.nitr            = v.taunit.*10;             % yr; nitrification timescale TUNED!
    tau.oxrm.s          = 0.05.*v.tauoxrm;          % yr; oxic remin rate (r = 0.0302 /day from Kriest + Oschlies, 2008) - tuned
    tau.oxrm.d          = 5.*v.tauoxrm;             % yr; (as above)
    tau.oxrm.n          = 2.*tau.sink.n.*v.tauoxrm; % yr; OM respiration rates vary from weeks to years, depending on rain rate (Martin + Sayles, 2004) - assume high rain rate but lower than open ocean
    tau.oxrm.z          = 2.*tau.sink.z.*v.tauoxrm; % yr; OM respiration rates vary from weeks to years, depending on rain rate (Martin + Sayles, 2004) - assume 100x lower rain rate than shallow seds
    tau.diss.s          = 0.25;                     % yr; CaCO3 dissolution timescale (tuned)
    tau.diss.d          = 0.25;                     % yr; CaCO3 dissolution timescale (tuned)
    tau.diss.n          = 10;                       % yr; CaCO3 dissolution timescale (tuned)
    tau.diss.z          = 100;                      % yr; CaCO3 dissolution timescale (tuned)
    tau.precip          = 25;           	        % yr; CaCO3 precipitation timescale (tuned)
    tau.photox          = 8.*(1/60).*(1/24).*(1/365);% yr; reduced iron abiotic photo-oxidation timescale (based on half-reaction rate at pH = 8 of 2-3 mins, Millero 1987) TUNED!
    tau.sorb            = 1e-1;                     % yr; slow reaction kinetics of phosphate iron sorption at T ~ 4C (Jaisi et al. 2010 experiment at this temp found 80% sorption over 1000 hrs, see fig 4)
    tau.mtrophy         = 1e-1;                     % yr; methanotrophy timescale, tuned to achieve modern CH4 air-sea flux (~10e12 g/yr from oceans, Cicero + Oremland, 1988 == 5e11 mol/yr)
    tau.mantle          = 1e9;                      % yr; mantle overturning timescale 
    tau.subd            = inp.sub;                  % yr; modern subduction timescale == lifetime of continental crust

end
