%% ============================= EONS Model =============================== 
% Julia Horne, 2018

% This script contains an assortment of literature references used to check
% the model's tuning. Everything is essentially for the modern Earth
% system, so one should not use these values to assess model output before 
% 4.5e9 yr
ConstantsParameters;

%% Organizational notes
% literature values for fluxes and reservoirs are collected below. The
% structure for these values denotes several important factors:
% ------- Reservoir naming : (rr,rp).RESERVOIR.SPECIES.AUTHOR.UNITS -------
% ----- Flux naming : (fr,fp).FLUXNAME.RESERVOIR.SPECIES.AUTHOR.UNITS -----
% rr,fr = reservoir/flux range, rp,fp = reservoir/flux preferred value 
% (any preferred values without noted error ranges are given a 25% error range)

% UNITS can take the forms:
% umolar  = µM           == micro molar     = 1e-6 mol/kg
% gcm2yr  = gC/m^2/yr    == grams C/ m^2/yr
% pgc     = PgC          == petagrams C     = 1e12 kg C
% tmol    = Tmol         == teramole        = 1e12 mol C
% pgcyr   = PgC/yr       == petagrams C/yr  = 1e12 kg C/yr
% molm2yr = mol/m^2/yr 
% tmolyr  = Tmol/yr      == teramole /yr    = 1e12 mol /yr
% gtcyr   = Gt C/yr      == gigatonnes C/yr = 1e15 g C/yr 
% gccm2yr = gC/cm^2/yr   == grams C/cm^2/yr  
% molyr   = mol/yr     
% gcyr    = g C/yr       == grams C/yr
% tgn     = Tg N         == teragrams N
% tgnyr   = Tg N/yr      == teragrams N/yr
% ppb     = parts per bil== 1e-9 
% molkg   = mol/kg       

% RESERVORS can take the forms:
% s = surface ocean
% d = deep ocean
% z = deep reactive sediments
% n = shelf reactive sediments
% u = unreactive sediments
% c = continent
% a = atmosphere
% t = total (usually total ocean)


%% ------------------------------ CARBON ----------------------------------
% from Chemical Oceanography, The Carbonate System (Millero 1996)
rr.s.CO2.m96.umolar              = [2000, 2400];     % µM; surface TCO2
rr.d.CO2.m96.umolar              = [2180, 2450];     % µM; deep TCO2
rr.s.TA.m96.umolar               = [2300, 2425];     % µM; surface Carbonate Alk
rr.d.TA.m96.umolar               = [2400, 2500];     % µM; deep Carbonate Alk
rr.t.TA.m96.umolar               = rr.s.TA.m96.umolar + rr.d.TA.m96.umolar;  

% from Bauer et al. 1995
rr.z.DIC.b95.umolar              = [2329, 4012];     % µM; sedimentary DIC
fr.diff.z.DIC.b95.gcm2yr         = [4, 22];          % gC/m2/yr; diffusive flux of DIC out of sediments
fp.diff.z.DIC.b95.gcm2yr         = 13;               % gC/m2/yr; diffusive flux of DIC out of sediments

% Libes (Introduction to Marine Biogeochemistry, 1992)
rr.n.OC.l92.pgc                  = [11, 14];         % mass; sedimentary POC
rr.n.CaCO3.l92.pgc               = [49, 72];         %? mass ; sedimentary CaCO3
fr.sed.n.OC.l92.molm2yr          = [0.027, 0.05];    % mol/m2/yr; marine sedimentation rate POC
fr.sed.n.CaCO3.l92.molm2yr       = [0.042, 0.08];    % mol/m2/yr; marine sedimentation rate CaCO3

% Watson + Orr 2003 and Doney et al. 2003 (Ocean Biogeochemistry, CH 5 and Ch 9)
rr.a.CO2.wod03.pgc               = [590, 750];       % Pg C; atmospheric CO2
rr.s.LB.wod03.pgc                = [2, 3];           % Pg C; total biota (live)
rr.d.DIC.wod03.pgc               = [38000, 38100];   % Pg C; total deep ocean DIC
rr.s.DIC.wod03.pgc               = [850, 1020];      % Pg C; total surface ocean DIC
rp.a.CO2.wod03.pgc               = 590;              % Pg C; atmospheric CO2
rp.s.LB.wod03.pgc                = 2;                % Pg C; total biota (live)
rp.d.DIC.wod03.pgc               = 38000;            % Pg C; total deep ocean DIC
rp.s.DIC.wod03.pgc               = 850;              % Pg C; total surface ocean DIC
rp.z.OC.wod03.pgc                = 150;              % Pg C; total sed carbon
fp.gasex.a.CO2.wod03.pgcyr       = 50;               % Pg C/yr; gas exchange 
fp.wthr.c.carb.wod03.pgcyr       = 0.7;              % Pg C/yr; riverine carbonate transport
fp.export.s.POC.wod03.pgcyr      = 2;                % Pg C/yr; POC export
fp.sed.z.POC.wod03.pgcyr         = 0.2;              % Pg C/yr; POC sedimentation
fp.meta.c.CaCO3.wod03.pgcyr      = 0.2;              % Pg C/yr; CO2 outgassing from continent
fr.gasex.a.CO2.wod03.pgcyr       = [50, 92];         % Pg C/yr; gas exchange 
fr.wthr.c.CaCO3.wod03.pgcyr      = [0.7, 0.8];       % Pg C/yr; riverine carbonate transport
fr.export.s.POC.wod03.pgcyr      = [2, 8];           % Pg C/yr; POC export
fr.sed.z.POC.wod03.pgcyr         = [0.15, 0.25];     % Pg C/yr; POC sedimentation
fr.meta.c.CaCO3.wod03.pgcyr      = [0.15, 0.25];     % Pg C/yr; CO2 degassing from continents

% from Fundamentals of Geobiology ch 2(Falkowski, 2012)
rp.a.CO2.f12.gtc                 = 835;              % Gt C; atm CO2 
rp.s.DIC.f12.gtc                 = 670;              % Gt C; surf DIC
rp.d.DIC.f12.gtc                 = 36730;            % Gt C; deep DIC
rp.t.OC.f12.gtc                  = 1e3;              % Gt C; total ocean biology (LB + DB)
rp.t.DIC.f12.gtc                 = 37400;            % Gt C; total inorganic carbon
rp.c.CaCO3.f12.gtc               = 6e7;              % Gt C; sedimentary carbonate on continent
rp.c.OC.f12.gtc                  = 15e6;             % Gt C; organic carbon (kerogens)
rp.t.CaCO3.f12.gtc               = 2500;             % Gt C; total inorganic ocean C
rp.s.LB.f12.gtc                  = 1;                % Gt C; active photosynthetic biomass  
rr.s.LB.f12.gtc                  = [1, 2];           % Gt C; aquatic biosphere
fr.netcp.t.CO2.f12.gtcyr         = [40, 50];         % Gt C/yr; total net carbon fixation by ocean biomass
fp.netcp.t.CO2.f12.gtcyr         = 45;               % Gt C/yr; total net carbon fixation by ocean biomass

% from Encyclopedia of Geochemistry (Carbon Cycle, Canuel and Hardison,2018)
rp.s.DIC.ch18.pgc                = 900;            	 % Pg C; surf DIC
rp.d.DIC.ch18.pgc                = 37100;            % Pg C; deep DIC
rp.u.CaCO3.ch18.pgc              = 1750;             % Pg C; sedimentary carbonate
rp.s.LB.ch18.pgc                 = 3;                % Pg C; total ocean biology (LB)
rp.a.CO2.ch18.pgc                = 589;              % Pg C; total atm CO2 
fr.gasex.a.CO2.ch18.pgcyr        = [60, 60.7];       % Pg C/yr; air-sea gas exchange CO2
fr.prod.s.CO2.ch18.pgcyr         = [37, 50];         % Pg C/yr; biological CO2 fixing
fr.ammon.d.CO2.ch18.pgcyr        = [2, 11];      	 % Pg C/yr; deep ocean CO2 remineralization
fp.gasex.a.CO2.ch18.pgcyr        = 60;           	 % Pg C/yr; air-sea gas exchange CO2
fp.sink.d.OC.ch18.pgcyr          = 90;           	 % Pg C/yr; DIC sinking to deep ocean
fp.burial.z.OC.ch18.pgcyr        = 0.2;              % Pg C/yr; POC burial
fp.prod.s.CO2.ch18.pgcyr         = 50;               % Pg C/yr; CO2 fixing
fp.ammon.s.CO2.ch18.pgcyr        = 37;               % Pg C/yr; surf ocean CO2 remineralization
fp.ammon.d.CO2.ch18.pgcyr        = 11;               % Pg C/yr; deep ocean CO2 remineralization
fp.ammon.t.CO2.ch18.pgcyr        = 48;               % Pg C/yr; total CO2 remineralization (s+d)
fp.wthr.c.sil.ch18.pgcyr         = 0.9;              % Pg C/yr; riverine HCO3 transport
fp.volc.u.CO2.ch18.pgcyr         = 0.1;          	 % Pg C/yr; volcanic CO2 outgassing-- assume total inorganic + organic

% from IPCC Ar5 WG1 Ch 6 (Carbon and other Biogeochemical Cycles, 2014)
rr.a.CH4.ar14.ppb                = [600, 700];       % ppb; modern atm CH4
rp.a.CH4.ar14.ppb                = 700;              % ppb; modern atm CH4 
rp.t.OC.ar14.pgc                 = 700;              % Pg C; dissolved organic carbon (LB+DB)
rp.t.DIC.ar14.pgc                = 38e3;             % Pg C; total ocean DIC 
% tau.POC.ar14.yr                  = 1e3;              % yr; minimum turnover time for POC
% tau.death.ar14.yr                = 2/52;             % yr; maximum turnover time for LB (2 weeks)
% tau.diss.ar14.yr                 = 10e3;             % yr; minimum timescale of CaCO3 dissolution --> 2HCO3
% tau.wthr.sil.ar14.pgc            = 10e6;             % yr; maximum timescale of CaSiO3 weathering
fp.gasex.a.CO2.ar14.pgcyr        = 80;               % Pg C/yr; air-sea gas exchange CO2
fp.volc.u.CO2.ar14.pgcyr         = 0.1;              % Pg C/yr; volcanic CO2 outgassing
fp.export.s.OC.ar14.pgcyr        = 90;               % Pg C/yr; productivity export
fp.burial.n.CaCO3.ar14.pgcyr     = 0.2;              % Pg C/yr; DIC burial
fp.remin.s.CO2.ar14.pgcyr        = 37;               % Pg C/yr; surface remin
fp.prod.s.CO2.ar14.pgcyr         = 50;               % Pg C/yr; CO2 assimilation 
fp.death.s.LB.ar14.pgcyr         = 2;                % Pg C/yr; death flux
fp.remin.d.CO2.ar14.pgcyr        = 13;               % Pg C/yr; deep remin
fp.remin.t.CO2.ar14.pgcyr        = 50;               % Pg C/yr; total remineralization (s+d)
fp.methox.t.CH4.ar14.molyr       = 3.2e13;           % mol C/yr: methane oxidation flux
fr.methox.t.CH4.ar14.molyr       = [0.75*3.2, 1.25*3.2].*1e13; % mol C/yr; methane oxidation flux
fp.netcp.t.CO2.ar14.pgcyr        = 90 - 50;          % Pg C/yr; net primary productivity, calculated (prod export - total remineralization)

% from Fundamentals of Geobiology ch 3(Wallmann + Aloisi, 2012 table 3.4 == balanced fluxes and references therein)
rp.c.CaCO3.wa12.mol              = 5e21;             % mol C; carbonate in continents
rp.c.OC.wa12.mol                 = 1.25e21;          % mol C; organic carbon in continents
% rr.z.CH4.wa12.mol                = [0.1, 0.3].*1e18; % mol C; CH4 in methane hydrates (assumed on seafloor) % this doesn't work because we don't resolve this kind of sequestration mechanism
rp.d.DIC.wa12.mol                = 2.7e18;           % mol C; deep ocean DIC
rp.s.DIC.wa12.mol                = 1e17;             % mol C; surface ocean DIC
rp.a.CO2.wa12.mol                = 6e15;             % mol C; atmospheric CO2 pre-industrial
rp.s.LB.wa12.mol                 = 5e14;             % mol C; total ocean biosphere (LB ?)
fr.mantle.t.CO2.wa12.tmolyr      = [3.1, 5.5];       % Tmol C/yr; mantle carbon outgassing (they assume CO2... spreading center release)
fp.mantle.t.CO2.wa12.tmolyr      = 4.3;              % Tmol C/yr; mantle carbon outgassing
fr.meta.c.CaCO3.wa12.tmolyr      = [2, 4];           % Tmol C/yr; CO2 flux from carbonate metamorphism
fp.meta.c.CaCO3.wa12.tmolyr      = 2.5;              % Tmol C/yr; CO2 flux from carbonate metamorphism
fr.wthr.c.CaCO3.wa12.tmolyr      = [10, 16];         % Tmol C/yr; carbonate chemical weathering
fp.wthr.c.CaCO3.wa12.tmolyr      = 11.7;          	 % Tmol C/yr; carbonate chemical weathering 
fr.meta.c.OC.wa12.tmolyr         = [0.4, 0.6];       % Tmol C/yr; CO2 flux from organic carbon metamorphism
fp.meta.c.OC.wa12.tmolyr         = 0.5;              % Tmol C/yr; CO2 flux from organic C metamorphism
fr.wthr.c.OC.wa12.tmolyr         = [8, 16];          % Tmol C/yr; weathering of continental organic matter
fp.wthr.c.OC.wa12.tmolyr         = 9;                % Tmol C/yr; weathering of continental organic matter
fr.burial.z.CaCO3.wa12.tmolyr    = [14, 17];         % Tmol C/yr; burial of biogenic carbonate at seafloor
fp.burial.z.CaCO3.wa12.tmolyr    = 16;               % Tmol C/yr; burial of biogenic carbonate at seafloor
fr.burial.z.OC.wa12.tmolyr       = [5.4, 27];      	 % Tmol C/yr; organic carbon burial in deep ocean
fp.burial.z.OC.wa12.tmolyr       = 10;               % Tmol C/yr; POC burial
fr.wthr.c.sil.wa12.tmolyr        = [6, 10];       	 % Tmol C/yr; weathering of silicate rocks
fp.wthr.c.sil.wa12.tmolyr        = 7.1;              % Tmol C/yr; weathering of silicate rocks
fp.export.s.POC.wa12.tmolyr      = 800;              % Tmol C/yr; export of POC to deep ocean
fr.volc.u.CO2.wa12.tmolyr        = [0.3, 0.5];       % Tmol C/yr; volcanic CO2 outgassing (total from mantle == 3.1 - 5.5 Tmol/yr)
fr.export.s.CaCO3.wa12.tmolyr    = [40, 130];        % Tmol C/yr; export of calcite
fr.diss.t.CaCO3.wa12.tmolyr      = [30, 120];        % Tmol C/yr; carbonate dissolution in deep ocean (export - 10 Tmol)
fp.prod.s.LB.wa12.gtcyr          = 48.5;             % Gt C/yr; net primary productivity
fp.prod.s.CO2.wa12.tmolyr        = 4040;             % Tmol C/yr; CO2 fixation 
fp.sink.d.OC.wa12.tmolyr         = 190;              % Tmol C/yr; organic C rain onto seafloor (ie. sinking velocity)
fp.volc.u.CO2.wa12.tmolyr        = 0.4;          	 % Tmol C/yr; volcanic CO2 outgassing (total from mantle ~ 4.3 Tmol/yr)
fp.precip.t.CaCO3.wa12.tmolyr    = 85;               % Tmol C/yr; carbonate precipitation in open oceans
fp.diss.t.CaCO3.wa12.tmolyr      = 80;               % Tmol C/yr; carbonate deep dissolution (majority of exported CaCO3)

% Wallman et al. 2008
fr.revweather.t.TA.w08.tmolyr    = [7.6, 28.8];      % Tmol C/yr; reverse weathering (Alk sink)
fr.sfw.t.CaCO3.w08.tmolyr        = [5, 20];          % Tmol C/yr; seafloor silicate weathering (DIC sink) 
fp.revweather.t.TA.w08.tmolyr    = 18.2;             % Tmol C/yr; reverse weathering (alk sink)

% Stumm + Morgan 1981 (Aquatic Chemistry)
fr.sed.z.OC.sm81.gccm2yr         =  [6e-4, 6e-1];    % g/cm2/yr; sedimentation rate of particles 
fr.diff.z.CO2.sm81.gccm2yr       =  [1e-5, 1e-1];    % g/cm2/yr; diffusion across sedimentary interface 

% Li + Elderfield, 2013 and references therein (5% error est)
fr.wthr.c.CaCO3.le13.molyr       = [11.7 12.9].*1e12;% mol C/yr; carbonate weathering
fp.wthr.c.CaCO3.le13.molyr       = 12.3e12;          % mol C/yr; carbonate weathering
fp.volc.u.CO2.le13.molyr         = 6e12;             % mol C/yr; CO2 outgassing from volcansim -- assume total inorganic + organic
fr.volc.u.CO2.le13.molyr         = [5.7 6.3].*1e12;  % mol C/yr; CO2 outgassing from volcansim
fr.wthr.c.sil.le13.molyr         = [8.27 9.13].*1e12;% mol C/yr; HCO3 produced in silicate weathering
fp.wthr.c.sil.le13.molyr         = 8.7e12;           % mol C/yr; HCO3 produced in silicate weathering
fp.burial.n.CaCO3.le13.molyr     = 16.8e12;          % mol C/yr; carbonate burial
fp.burial.z.CaCO3.le13.molyr     = 16.8e11;          % mol C/yr; carbonate burial (assuming most is on cont shelf)
fp.wthr.c.OC.le13.molyr          = 3.6e12;           % mol C/yr; oxidative weathering of Org C
fp.burial.z.OC.le13.molyr        = 5.1e12;           % mol C/yr; organic C burial
fp.revweather.t.TA.le13.molyr    = 4.6e12;           % mol C/yr; HCO3 consumed by reverse weathering
fr.revweather.t.TA.le13.molyr    = [4.37 4.83].*1e12;% mol C/yr; HCO3 consumed by reverse weathering

% from Rigwell et al 2007
fp.export.s.CaCO3.r07.pgcyr      = 1.2;              % Pg C/yr; export production of CaCO3
fp.export.s.OC.r07.pgcyr         = 8.9;              % Pg C/yr; export production of DB

% Bergman et al 2004 COPSE
rr.c.CaCO3.b04.mol               = [9.375e20, 1.5625e21];% mol C; total continental inorganic carbon
rr.c.OC.b04.mol                  = [3.75e21, 6.25e21]; % mol C; total continental organic carbon
rp.c.CaCO3.b04.mol               = 1.25e21;          % mol C; total continental inorganic carbon
rp.c.OC.b04.mol                  = 5e21;             % mol C; total continental organic carbon
fp.wthr.c.OC.b04.molyr           = 7.75e12;          % mol C/yr; oxidative erosion
fp.meta.c.OC.b04.molyr           = 1.125e12;         % mol C/yr; organic carbon degassing
fp.mantle.t.CH4.b04.molyr        = 7.5e10;           % mol C/yr; reduced degassing flux from the mantle (Goldblatt et al. 2006 - 7.5e10 mol O2/yr in modern)
fr.mantle.t.CH4.b04.molyr        = [1.25, 7.5].*1e10;% mol C/yr; mantle degassing flux

% % Holland 1978 via Pavlov et al. 2001 
% fp.haze.a.HZ.h78.gcyr            = 1e14;             % g C/yr; photochemical haze production at pCO2 ~ 2500 ppm 

% Williams et al. 1992 
fp.volc.u.CO2.w92.gcyr           = 65e12;            % g C/yr; global CO2 volcanic flux

% Zeebe 2012 LOSCAR and references therein (Walker + Kasting, 1992; Morse + Mackenzie, 1990)
fp.wthr.c.sil.z12.molyr          = 5e12;             % mol C/yr; silicate weathering flux, modern 
fp.wthr.c.CaCO3.z12.molyr        = 12e12;            % mol C/yr; carbonate weathering flux, modern

% Feely et al. 2004 + references therein (Feely 2002)
fp.remin.s.OC.f04.pgcyr          = 5.3;              % Pg C/yr; upper water column remineralization flux
fr.remin.s.OC.f04.pgcyr          = [4.3, 6.3];       % Pg C/yr; upper water column remineralization flux
fp.diss.t.CaCO3.f04.pgcyr        = 0.31;             % Pg C/yr; total ocean CaCO3 remineralization flux (dissolution?)

% Martin + Sayles, 2004
fr.ammon.n.OC.ms04.molyr         = [50, 56].*1e5.*v.sed.sa.n.*1e-6; % mol/yr; oxygen consumption rate converted from  umol/cm2/yr
fp.ammon.n.OC.ms04.molyr         = 53.*1e5.*v.sed.sa.n.*1e-6; % mol/yr; oxygen consumption rate converted from  umol/cm2/yr

% Wong et al. 2019 + references therein (converted from tonnes/megatonnes to gigatonnes) 
fr.subduct.u.CaCO3.w19.gtcyr     = [13 57].*1e-3;   % Gt C/yr; carbonate sediment subduction flux (mix of Kelemen + Manning 2015 and Dutkewicz et al. 2018 for low + high ests)
fr.subduct.o.CaCO3.w19.gtcyr     = [20 55].*1e-3;   % Gt C/yr; carbonate sediment subduction flux (fig. 1)
fr.volc.u.CO2.w19.gtcyr          = [70 100].*1e-3;  % Gt C/yr; CO2 volcanic flux estimate range -- assume total inorganic + organic
fr.volc.u.carb.w19.gtcyr         = fr.volc.u.CO2.w19.gtcyr.*0.8; % Gt C/yr; assuming that 80% downgoing carbon is carbonate
fr.volc.u.OC.w19.gtcyr           = fr.volc.u.CO2.w19.gtcyr.*0.2; % Gt C/yr; assuming that 80% downgoing carbon is carbonate
fp.mantle.t.CO2.w19.gtcyr        = 21.*1e-3;        % Gt C/yr; estimate for MORB CO2 degassing (plume volc is 10x lower, and cannot constrain C source from mantle/crust, so MORB is taken as mantle source)
fr.mantle.t.CO2.w19.gtcyr        = [8 33].*1e-3;    % Gt C/yr; ^^ 21±13
rr.u.CaCO3.w19.gtc               = fr.subduct.u.CaCO3.w19.gtcyr.*1e8;% Gt C; assuming ALL sediments are subducted, the higher-end estimate for pelagic sedimentary carbonate
rr.o.CaCO3.w19.gtc               = fr.subduct.o.CaCO3.w19.gtcyr.*1e8;% Gt C; assuming ALL crust is subducted, the higher-end estimate for oceanic crust carbonate
rp.u.CaCO3.w19.gtc               = mean(rr.u.CaCO3.w19.gtc); % for want of a prefered reference...
rr.m.C.w19.pgc                   = [636,726].*1e-6.*v.m.mass./1e12;% Pg C; serpentinized mantle estimate in ppm converted to Pg C - Kelemen + Manning, 2015
rp.m.C.w19.pgc                   = 681.*1e-6.*v.m.mass./1e12;% Pg C; serpentinized mantle estimate in ppm converted to Pg C

% Clift 2017
fp.subduct.u.C.c17.gtcyr         = 60.*1e-3;        % Gt C/yr; total sediment carbon subduction estimate
fr.subduct.u.C.c17.gtcyr         = [48 72].*1e-3;   % Gt C/yr; with 20% error estimate
rp.u.OC.c17.gtc                  = fp.subduct.u.C.c17.gtcyr.*1e8.*0.2; % Gt C; 20 % total downgoing carbon assumed to be organic 
rr.u.OC.c17.gtc                  = fr.subduct.u.C.c17.gtcyr.*1e8.*0.2; % Gt C; 20 % total downgoing carbon assumed to be organic 
rp.u.CaCO3.c17.gtc               = fp.subduct.u.C.c17.gtcyr.*1e8.*0.8; % Gt C; 80 % total downgoing carbon assumed to be carbonate 
rr.u.CaCO3.c17.gtc               = fr.subduct.u.C.c17.gtcyr.*1e8.*0.8; % Gt C; 80 % total downgoing carbon assumed to be carbonate 

% Mackenzie et al. 2004 + references therein (table 1 and figs 1+2) - converted from grams/tons C to pgC (/1e15, /1e9) 
rr.m.C.m04.pgc                   = [8.9 16.6].*1e7;  % Pg C; upper mantle estimate taken from Li 2000 and Wood 1996
rp.o.CaCO3.m04.pgc               = 9.2e5;            % Pg C; oceanic crust carbon (assumed all carbonate) 
rp.t.DIC.m04.pgc                 = 3.85e4;           % Pg C; total ocean carbon 
rp.t.CaCO3.m04.pgc               = 6.53e7;           % Pg C; total carbonates - assumed to be total ocean carbonates (?)
rp.n.CaCO3.m04.pgc               = 6.53e7/2;         % Pg C; shelf carbonates - assuming 50% of buried ocean carbonate is on shelf (pg. 20)
rp.z.CaCO3.m04.pgc               = 6.53e7/2;         % Pg C; seafloor carbonates - see above
fp.sed.z.CaCO3.m04.pgcyr         = 3.2;              % Pg C; carbonate sedimentation rate at seafloor
fp.sed.z.OC.m04.pgcyr            = 0.5;              % Pg C; organic sedimentation rate at seafloor
fr.sed.n.CaCO3.m04.molyr         = [8 12.5].*1e12;   % mol/yr; sedimentation rate of inorg C at shelf
fr.sed.z.CaCO3.m04.molyr         = [8.8 11.6].*1e12; % mol/yr; sedimentation rate of inorg C at slope and open ocean
fr.precip.n.CaCO3.m04.molyr      = [9 15.5].*1e12;   % mol/yr; precipitated CaCO3 on shelf
fr.precip.s.CaCO3.m04.molyr      = [4.8 65].*1e12;   % mol/yr; precipitated CaCO3 in slope/ocean
fr.precip.t.CaCO3.m04.molyr      = fr.precip.s.CaCO3.m04.molyr + fr.precip.n.CaCO3.m04.molyr; 
fr.diss.d.CaCO3.m04.molyr        = [8.8 54].*1e12;   % mol/yr; dissolved CaCO3 in slope/ocean
fr.burial.z.CaCO3.m04.molyr      = [6 11.6].*1e12;   % mol/yr; slope + seafloor accumulated (buried?) CaCO3
fr.burial.n.CaCO3.m04.molyr      = [7 7.5].*1e12;    % mol/yr; shelf accumulated (buried?) CaCO3
fr.ammon.n.OC.m04.molyr          = [30.4 93].*1e12;  % mol/yr; remin (anox + ox) OC on shelf proximal + distal - assuming most is oxic!
fr.ammon.z.OC.m04.molyr          = [34.7 65].*1e12;  % mol/yr; remin (anox + ox) OC on slope and open ocean - assuming most is oxic!
fr.sed.n.OC.m04.molyr            = [29.4 100].*1e12; % mol/yr; sedimentation OC on shelf proximal + distal - assuming most is oxic!
fr.sed.z.OC.m04.molyr            = [35 67].*1e12;    % mol/yr; sedimentation OC on slope and open ocean - assuming most is oxic!
fr.burial.n.OC.m04.molyr         = [6.5 7.5].*1e12;  % mol/yr; burial OC on shelf proximal + distal - assuming most is oxic!
fr.burial.z.OC.m04.molyr         = [0.3 2].*1e12;    % mol/yr; burial OC on slope and open ocean - assuming most is oxic!
rr.n.OC.m04.mol                  = fr.burial.n.OC.m04.molyr.*1e2;% mol; reservoir of OC in shallow seds, assuming a flux mutiplied by the average burial timescale in modern
rr.z.OC.m04.mol                  = fr.burial.z.OC.m04.molyr.*1e3;% mol; reservoir of OC in deep seds, assuming a flux mutiplied by the average burial timescale in modern
fr.diss.t.CaCO3.m04.molyr        = fr.precip.t.CaCO3.m04.molyr - (fr.burial.z.CaCO3.m04.molyr + fr.burial.n.CaCO3.m04.molyr); % total dissolution assumed to be difference in precip and burial


% all values without author-quoted ranges are given 25% error bars
rr.a.CO2.f12.gtc                 = [626, 1045];      % Gt C; atm CO2
rr.s.DIC.f12.gtc                 = [503, 837.5];     % Gt C; surf DIC
rr.d.DIC.f12.gtc                 = [2.755e4,4.6e4];  % Gt C; deep DIC
rr.t.OC.f12.gtc                  = [750, 1250];      % Gt C; total ocean biology (LB + DB)
rr.c.CaCO3.f12.gtc               = [4.5e7, 10.5e7];  % Gt C; continental carbonate
rr.c.OC.f12.gtc                  = [11.25e6, 18.75e6];% Gt C; continental organic carbon
rr.t.CaCO3.f12.gtc               = [625 3125];       % Gt C; total ocean CaCO3
rr.s.DIC.ch18.pgc                = [675, 1125];   	 % Pg C; surf DIC
rr.d.DIC.ch18.pgc                = [27825, 46375]; 	 % Pg C; deep DIC
rr.u.CaCO3.ch18.pgc              = [1.3e3,2.2e3];    % Pg C; continental carbonate
rr.s.LB.ch18.pgc                 = [2.25, 3.75];     % Pg C; total ocean biology (LB)
rr.a.CO2.ch18.pgc                = [442, 736.25];    % Pg C; total atm CO2 
rr.t.OC.ar14.pgc                 = [575, 875];       % Pg C; total ocean biology (LB + DB)
rr.t.DIC.ar14.pgc                = [28.5e3, 47.5e3]; % Pg C; total ocean DIC  
rr.u.OC.wod03.pgc                = [112.5, 187.5];   % Pg C; total sed carbon
rr.s.LB.wa12.mol                 = [3.75, 6.25].*1e14;% mol C; surface ocean LB 
fr.prod.s.LB.wa12.gtcyr          = [36.4, 60.6];     % Gt C/yr; net primary productivity
fr.prod.s.CO2.wa12.tmolyr        = [3030, 5050];     % Tmol C/yr; CO2 fixation 
fr.precip.n.CaCO3.wa12.tmolyr    = [7.5, 12.5];      % Tmol C/yr; carbonate precipitation in shelf ocean
fr.precip.t.CaCO3.wa12.tmolyr    = [40, 130];        % Tmol C/yr; carbonate precipitation in open ocean
fr.prod.s.CO2.ch18.pgcyr         = [12.5, 62.5];     % Pg C/yr; CO2 assim
fr.sink.d.OC.ch18.pgcyr          = [67.5, 112.5];  	 % Pg C/yr; OC sinking to deep ocean
fr.burial.n.CaCO3.ch18.pgcyr     = [0.15, 0.25];     % Pg C/yr; DIC burial
fr.ammon.s.CO2.ch18.pgcyr        = [27.8, 46.25];    % Pg C/yr; surf ocean CO2 remineralization
fr.ammon.t.CO2.ch18.pgcyr        = [36, 60];         % Pg C/yr; total CO2 remineralization (s+d)
fr.wthr.c.sil.ch18.pgcyr         = [0.675, 1.125];   % Pg C/yr; riverine HCO3 transport
fr.volc.u.CO2.ch18.pgcyr         = [0.075, 0.125];   % Pg C/yr; volcanic CO2 outgassing -- assume total inorganic + organic
fr.netcp.t.CO2.ch18.pgcyr        = [1.875, 19.875];  % Pg C/yr; calculated net primary productivity (highest est. total CO2 prod - total remin ) 
fr.gasex.a.CO2.ar14.pgcyr        = [60, 140];        % Pg C/yr; CO2 gas exchange
fr.volc.u.CO2.ar14.pgcyr         = [0.075, 0.125];   % Pg C/yr; CO2 volcanic outgassing
fr.export.s.OC.ar14.pgcyr        = [67.5, 112.5];    % Pg C/yr; export production
fr.burial.n.CaCO3.ar14.pgcyr     = [0.15, 0.25];     % Pg C/yr; DIC burial
fr.remin.s.CO2.ar14.pgcyr        = [27.75, 46.25];   % Pg C/yr; surf remin
fr.remin.d.CO2.ar14.pgcyr        = [9.75, 16.25];    % Pg C/yr; deep remin
fr.remin.t.CO2.ar14.pgcyr        = [37.5, 62.5];     % Pg C/yr; total remin
fr.prod.s.CO2.ar14.pgcyr         = [37.5, 62.5];     % Pg C/yr; CO2 assimilation
fr.death.s.LB.ar14.pgcyr         = [1.5, 2.5];       % Pg C/yr; death flux
fr.netcp.t.CO2.ar14.pgcyr        = [30, 50];         % Pg C/yr; net primary productivity, calculated (prod export - total remineralization)
fr.export.s.CaCO3.r07.pgcyr      = [0.9, 1.5];       % Pg C/yr; CaCO3 export 
fr.export.s.OC.r07.pgcyr         = [6.675, 11.125];  % Pg C/yr; POC export
fr.wthr.c.OC.b04.molyr           = [5.8e12, 9.68e12];% mol C/yr; oxidative weathering
fr.meta.c.OC.b04.molyr           = [0.85e12, 1.4e12];% mol C/yr; organic c degassing
fr.sed.z.OC.wod03.pgcyr          = [0.15, 0.25];     % Pg C/yr; POC sedimentation
fr.meta.c.CO2.wod03.pgcyr        = [0.15, 0.25];     % Pg C/yr; CO2 degassing from continents
% fr.haze.a.HZ.h78.gcyr            = [7.5, 12.5].*1e13;% g C/yr; photochemical haze production
fr.burial.n.CaCO3.le13.molyr     = [12.6, 21].*1e12; % mol C/yr; carbonate burial (assuming most in on cont shelf)
fr.burial.z.CaCO3.le13.molyr     = 0.1.*fr.burial.n.CaCO3.le13.molyr; % mol C/yr; carbonate burial ^^
fr.wthr.c.OC.le13.molyr          = [2.7, 4.5].*1e12; % mol C/yr; oxidative weathering of Org C
fr.burial.z.OC.le13.molyr        = [3.825, 6.375].*1e12;% mol C/yr; organic C burial
fr.wthr.c.sil.z12.molyr          = [3.75, 6.25].*1e12;% mol C/yr; silicate weathering flux
fr.wthr.c.CaCO3.z12.molyr        = [9, 15].*1e12;    % mol C/yr; carbonate weathering flux
fr.diss.t.CaCO3.f04.pgcyr        = [0.2325,0.3875];  % Pg C/yr; total ocean CaCO3 remineralization flux (dissolution?)
rr.s.LB.f12.gtc                  = [0.75, 1.25];     % Gt C; active photosynthetic biomass  
rr.o.CaCO3.m04.pgc               = [6.9 11.5].*1e5;  % Pg C; oceanic crust carbon (assumed all carbonate) 
rr.n.CaCO3.m04.pgc               = [3.06 4.08].*1e7; % Pg C; shelf carbonates - assuming 50% of buried ocean carbonate is on shelf (pg. 20)
rr.z.CaCO3.m04.pgc               = rr.n.CaCO3.m04.pgc;% Pg C; seafloor carbonates - see above


% assume that organic carbon reservoirs in continent and sediments reflect
% the organic N and P reservoirs!
rp.c.OP.f12.gtp                  = rp.c.OC.f12.gtc./v.const.CPratio; % Gt P; organic phosphate
rr.c.OP.f12.gtp                  = rr.c.OC.f12.gtc./v.const.CPratio; % Gt P; organic phosphate
rp.c.ON.f12.gtn                  = rp.c.OC.f12.gtc./v.const.CNratio; % Gt N; organic ammonia
rr.c.ON.f12.gtn                  = rr.c.OC.f12.gtc./v.const.CNratio; % Gt N; organic ammonia
rp.c.OP.wa12.mol                 = rp.c.OC.wa12.mol./v.const.CPratio;% mol P; organic phosphate
rp.c.ON.wa12.mol                 = rp.c.OC.wa12.mol./v.const.CNratio;% mol N; organic ammonia
rp.c.OP.b04.mol                  = rp.c.OC.b04.mol./v.const.CPratio; % mol P; organic phosphate
rr.c.OP.b04.mol                  = rr.c.OC.b04.mol./v.const.CPratio; % mol P; organic phosphate
rp.c.ON.b04.mol                  = rp.c.OC.b04.mol./v.const.CNratio; % mol N; organic ammonia
rr.c.ON.b04.mol                  = rr.c.OC.b04.mol./v.const.CNratio; % mol N; organic ammonia

%% ----------------------------- NITROGEN  --------------------------------

% from Fundamentals of Geobiology ch 4(Ward, 2012)
rp.a.N2.w12.tgn                  = 3.7e9;            % Tg N; atm N2
rp.t.N2.w12.tgn                  = 1.46e6;           % Tg N; total ocean N2
rp.t.HNO3.w12.tgn                = 6e5;              % Tg N; total ocean NO3
rp.t.PON.w12.tgn                 = 9e4;              % Tg N; total ocean biology (LB + DB)
rp.t.N.w12.tgn                   = 4.33e6;           % Tg N; total inorganic N
rp.c.NH4.w12.tgn                 = 1.3e9;            % Tg N; total continental crust N
fr.fix.s.N2.w12.tgnyr            = [110, 158];     	 % Tg N/yr; total oceanic N2 fixation
fr.denit.t.HNO3.w12.tgnyr        = [123, 2030];      % Tg N/yr; oceanic denitrification 
fp.volc.u.N2.w12.tgnyr           = 5;                % Tg N/yr; volcanic N2 outgassing
fp.burial.z.NH4.w12.tgnyr        = 25;               % Tg N/yr; burial of NH4 

% from Nitrogen in the Marine Environment (Gruber, 2008)
rr.a.N2.g08.tgn                  = [2.8e9, 5.2e9];	 % Tg N; atmospheric inventory of N2
rr.c.NH4.g08.tgn                 = [84e3, 156e3];	 % Tg N; continental inventory 
rr.t.PON.g08.tgn                 = [200, 600];       % Tg N; ocean particulate org. N (total biomass)
rr.t.N.g08.tgn                   = [7.4e6, 1.375e7]; % Tg N; total inorganic N 
rr.s.RN.g08.umolar               = [0.1, 1.2];       % µmol/kg; surface ocean NH4
rr.d.RN.g08.umolar               = [0.1, 0.55];      % µmol/kg; deep ocean NH4
rp.a.N2.g08.tgn                  = 4e9;            	 % Tg N; atmospheric inventory of N2
rp.c.NH4.g08.tgn                 = 1.2e5;         	 % Tg N; continental inventory 
rp.t.PON.g08.tgn                 = 400;            	 % Tg N; ocean particulate org. N (total biomass)
rp.t.N.g08.tgn                   = 10581250;         % Tg N; total inorganic N 
fr.volc.u.N2.g08.tgnyr           = [14,  26] ;       % Tg N/yr; volcanic N2 flux
fr.wthr.c.ON.g08.tgnyr           = [15, 80];  	     % Tg N/yr; river transport of NH4 (assumed mostly from Org matter)
fr.fix.s.N2.g08.tgnyr            = [70, 170] ;       % Tg N/yr; N2 fixation 
fr.export.s.ON.g08.tgnyr         = [350, 1120] ;     % Tg N/yr; export production 
fr.burial.z.NH4.g08.tgnyr        = [50, 350];        % Tg N/yr; burial flux 
fr.nitr.t.HNO3.g08.tgnyr         = [1025, 6200];     % Tg N/yr; nitrification flux
fr.denit.t.HNO3.g08.tgnyr        = [60, 220] ;       % Tg N/yr; denitrification flux 
fr.denit.s.HNO3.g08.tgnyr        = [5, 60];          % Tg N/yr; surf denit
fr.denit.d.HNO3.g08.tgnyr        = [5, 175];         % Tg N/yr; deep denit
fr.denit.z.HNO3.g08.tgnyr        = [5, 175];         % Tg N/yr; sedimentary denit
fr.assim.s.NH3.g08.tgnyr         = [1380, 7200];     % Tg N/yr; assimilation flux 
fr.ammon.s.NH3.g08.tgnyr         = [900, 6200];      % Tg N/yr; surf remin (either as NH4 or NO3, depending on O2 conc)
fr.ammon.d.NH3.g08.tgnyr         = [110, 900];       % Tg N/yr; deep remin
fr.ammon.t.NH3.g08.tgnyr         = [1025, 6200];     % Tg N/yr; remineralization flux
fp.volc.u.N2.g08.tgnyr           = 20;               % Tg N/yr; volcanic N2 flux
fp.wthr.c.ON.g08.tgnyr           = 15;               % Tg N/yr; river transport of NH4
fp.fix.s.N2.g08.tgnyr            = 120;              % Tg N/yr; N2 fixation 
fp.export.s.ON.g08.tgnyr         = 1120;             % Tg N/yr; export production 
fp.burial.z.NH4.g08.tgnyr        = 50;               % Tg N/yr; burial flux 
fp.nitr.t.HNO3.g08.tgnyr         = 1200;             % Tg N/yr; nitrification flux
fp.denit.t.HNO3.g08.tgnyr        = 85;               % Tg N/yr; denitrification flux 
fp.denit.s.HNO3.g08.tgnyr        = 65;               % Tg N/yr; surf denit
fp.denit.d.HNO3.g08.tgnyr        = 110;              % Tg N/yr; deep denit
fp.denit.z.HNO3.g08.tgnyr        = 5;                % Tg N/yr; sedimentary denit
fp.assim.s.NH3.g08.tgnyr         = 7000;             % Tg N/yr; assimilation flux 
fp.ammon.t.NH3.g08.tgnyr         = 6000;             % Tg N/yr; remineralization flux
fp.prod.s.NH3.g08.tgnyr          = 7120;             % Tg N/yr; fixing + assim = total N prod

% from WOCE Atlas online
rr.s.fixN.woce.umolar            = [5, 32];          % µmol N/kg ; total fixed N in surface ocean
rr.d.fixN.woce.umolar            = [30, 38];         % µmol N/kg ; total fixed N in deep ocean
rr.t.fixN.woce.umolar            = rr.s.fixN.woce.umolar + rr.d.fixN.woce.umolar;
rr.s.HNO3.woce.umolar            = [5, 42];          % µmol N/kg ; total nitrate in surface ocean
rr.d.HNO3.woce.umolar            = [28, 45];         % µmol N/kg ; total nitrate in deep ocean
rp.s.fixN.woce.umolar            = mean(rr.s.fixN.woce.umolar); 
rp.d.fixN.woce.umolar            = mean(rr.d.fixN.woce.umolar); 
rp.t.fixN.woce.umolar            = mean(rr.t.fixN.woce.umolar); 
rp.s.HNO3.woce.umolar            = mean(rr.s.HNO3.woce.umolar); 
rp.d.HNO3.woce.umolar            = mean(rr.d.HNO3.woce.umolar); 


% from Encyclopedia of Geochemistry (Nitrogen Cycle, Palta and Harnett, 2018)
fr.denit.t.HNO3.ph18.tgnyr       = [100, 280];    	 % Tg N/yr; oceanic denitrification
fp.fix.s.N2.ph18.tgnyr           = 140;              % Tg N/yr; N2 fixing
fp.gasex.a.NH3.ph18.tgnyr        = 9;                % Tg N/yr; NH3 ocean-atm emission
fp.burial.z.NH4.ph18.tgnyr       = 20;           	 % Tg N/yr; NH4 burial

% from Marine Nitrogen Cycle (Voss et al. 2013)
fr.wthr.c.ON.v13.tgnyr           = [40, 60];         % Tg N/yr; river NH4 transport (assumed mostly from Org Matter)
fr.denit.s.HNO3.v13.tgnyr        = [4, 100];         % Tg N/yr; surf denit
fr.denit.d.HNO3.v13.tgnyr        = [200, 300];       % Tg N/yr; deep denit
fr.fix.s.N2.v13.tgnyr            = [15, 177];        % Tg N/yr; N2 fixation
fp.prod.s.NH3.v13.tgnyr          = 320;              % Tg N/yr; net primary production
fp.fix.s.N2.v13.tgnyr            = 15;               % Tg N/yr; N2 fixation
fp.burial.z.NH4.v13.tgnyr        = 22;               % Tg N/yr; NH4 burial
fp.denit.s.HNO3.v13.tgnyr        = 60;               % Tg N/yr; surf denit
fp.denit.d.HNO3.v13.tgnyr        = 200;              % Tg N/yr; deep denit
fp.export.s.ON.v13.tgnyr         = 390;              % Tg N/yr; export production
fp.wthr.c.ON.v13.tgnyr           = 40;               % Tg N/yr; river NH4 transport

% from IPCC Ar5 WG1 Ch 6 (Carbon and other Biogeochemical Cycles, 2014)
fr.denit.t.HNO3.ar14.tgnyr       = [200, 400];       % Tg N/yr; total N2 release in denitrification
fr.fix.s.N2.ar14.tgnyr           = [140, 177];       % Tg N/yr; N2 fixation
fp.denit.t.HNO3.ar14.tgnyr       = 300;              % Tg N/yr; total N2 release in denitrification
fp.fix.s.N2.ar14.tgnyr           = 160;              % Tg N/yr; N2 fixation
fp.gasex.a.NH3.ar14.tgnyr        = 8.2;              % Tg N/yr; air-sea NH3 gas exchange 

% from COPSE model (Bergman et al. 2004)
fp.denit.t.HNO3.b04.molyr        = 8.6e12;           % mol N/yr; total denitrification
fp.fix.s.N2.b04.molyr            = 8.72e12;          % mol N/yr; total N2 fixation

% from Busigny et al. 2003 - Subduction zone recycling
fp.subduct.m.NH4.b03.molyr       = 7.6e11;           % mol N/yr; island arc subduction zone recycling estimate

% from Martin + Sayles 2004 (Treatise on Geochemistry, Vol. 7, ch. 2)
% note: burial rates avg. 1 cm/kyr and avg. particle remains in mixed
% (reactive) layer fo r~ 1e4 yr in deep ocean, ~1e2-1e3 yr in shallow seds
% fp.denit.z.NO3.ms04.tgnyr        = 230;             % Tg N/yr; deep sediment denitrification rates (over twice the rates of water column denit)

% from Goldblatt et al. (2009) - converted from 1e18 kg == 1e9 Tg 
fp.subduct.u.ON.g09.tgnyr        = 0.76;            % Tg N/yr; subduction of N in seds, assumed here to be organic material
fp.subduct.o.NH4.g09.tgnyr       = 0.12;            % Tg N/yr; subduction of N in altered and unaltered oceanic crust
fp.volc.u.ON.g09.tgnyr           = 0.55;            % Tg N/yr; volcanism of N, assumed here to be dominated by organic material from seds
rp.c.NH4.g09.tgn                 = 0.55e9;          % Tg N; continental reservoirs of igneous rock
rr.c.NH4.g09.tgn                 = [0.28 0.82].*1e9;% Tg N; continental reservoirs of igneous rock
rp.c.ON.g09.tgn                  = 1.55e9;          % Tg N; continental reservoirs of sed+metased, assumed mostly organic 
rr.c.ON.g09.tgn                  = [0.77 2.33].*1e9;% Tg N; continental reservoirs of sed+metased, assumed mostly organic
rp.o.NH4.g09.tgn                 = 12.4e6;          % Tg N; oceanic crust reservoirs (altered + not)
rr.o.NH4.g09.tgn                 = [6.2 18.6].*1e6; % Tg N; oceanic crust reservoirs (altered + not)
rp.u.ON.g09.tgn                  = 0.31e9;          % Tg N; oceanic crust sediments, assumed here to be mostly organic
rr.u.ON.g09.tgn                  = [0.15 0.47].*1e9; % Tg N; oceanic crust sediments, assumed here to be mostly organic
rp.m.N.g09.tgn                   = 8.4e9;           % Tg N; upper mantle reservoirs (MORB/OIB + sublithos)
rr.m.N.g09.tgn                   = [3.2 11.6].*1e9; % Tg N; upper mantle reservoirs (MORB/OIB + sublithos)

% from Byrne + Goldblatt, 2014 (table 1)
rp.a.NH3.bg14.mol                = 6.7e-10.*v.atm.mol;% mol N; atmospheric reservoir given modern mixing ratio ~ 7e-10 
fp.gasex.a.NH3.bg14.molyr        = 0.17.*1.3e12;     % mol N/yr; flux to atmosphere from sea surface (~ 17% of total flux)

% from Johnson + Goldblatt, 2015 (table 13) - preferred (bolded) values -- converted from 1e18 kg == 1e9 Tg 
rp.m.N.jg15.tgn                  = 24.2.*1e9;        % Tg N; total in mantle
rr.m.N.jg15.tgn                  = [8.2, 40.2].*1e9; % Tg N; total in mantle
rp.c.N.jg15.tgn                  = 2.8.*1e9;         % Tg N; total in BSE - mantle
rr.c.N.jg15.tgn                  = [1, 4.667].*1e9;  % Tg N; total in BSE - mantle

% all values without author-quoted ranges are given 25% error bars
rr.t.N.w12.tgn                   = [3.25,54.125].*1e5;% Tg N; total ocean inorganic N
rr.a.N2.w12.tgn                  = [2.8e9, 4.625e9]; % Tg N; atm N2
rr.t.N2.w12.tgn                  = [1.095e6, 1.8e6]; % Tg N; total ocean N2
rr.t.HNO3.w12.tgn                = [4.5e5, 7.5e5];   % Tg N; total ocean NO3
rr.c.NH4.w12.tgn                 = [9.75e8, 1.625e9];% Tg N; total continental crust N
rr.t.PON.w12.tgn                 = [675, 1125].*1e2; % Tg N; total organic N
fr.volc.u.N2.w12.tgnyr           = [3.75, 6.25];     % Tg N/yr; volcanic N2 outgassing
fr.burial.z.NH4.w12.tgnyr        = [18.75, 31.25];   % Tg N/yr; burial of NH4 
fr.fix.s.N2.ph18.tgnyr           = [105, 175];       % Tg N/yr; N2 fixing
fr.gasex.a.NH3.ph18.tgnyr        = [6.75, 11.25];	 % Tg N/yr; NH3 ocean-atm emission
fr.burial.z.NH4.ph18.tgnyr       = [15, 25];      	 % Tg N/yr; NH4 burial
fr.gasex.a.NH3.w12.tgnyr         = [6.15, 10.25];    % Tg N/yr; NH3 gas exchange
fr.burial.z.NH4.v13.tgnyr        = [16.5, 27.5];     % Tg N/yr; NH4 burial 
fr.prod.s.NH3.v13.tgnyr          = [240, 400];       % Tg N/yr; net primary production
fr.gasex.a.NH3.ar14.tgnyr        = [6.15, 10.25];    % Tg N/yr; NH3 gas exchange
fr.denit.z.NO3.ms04.tgnyr        = [172.5, 287.5];   % Tg N/yr; deep sediment denitrification rates (over twice the rates of water column denit)
fr.subduct.m.NH4.b03.molyr       = [5.7, 9.5] .*1e11;% mol N/yr; island arc subduction zone recycling estimate
fr.subduct.u.ON.g09.tgnyr        = [0.57 0.95];      % Tg N/yr; subduction of N in seds, assumed here to be organic material
fr.subduct.o.NH4.g09.tgnyr       = [0.09 0.15];      % Tg N/yr; subduction of N in altered and unaltered oceanic crust
fr.volc.u.ON.g09.tgnyr           = [0.412 0.687];    % Tg N/yr; volcanism of N, assumed here to be dominated by organic material from seds
rr.a.NH3.bg14.mol                = [5 8.4].*1e-10.*v.atm.mol;% mol N; atmospheric reservoir given modern mixing ratio ~ 7e-10 
fr.gasex.a.NH3.bg14.molyr        = [0.166 0.276].*1e12;% mol N/yr; flux to atmosphere from sea surface (~ 17% of total flux)


% from Nitrogen in the Marine Environment (Gruber, 2008)
RtauNO3                          = [350, 389];       % yr; turnover time
RtauNH4                          = [0.04, 0.06];     % yr; turnover time
RtauPON                          = [0.025, 0.075];   % yr; turnover time
RtauN2                           = [48600, 5.94e4];  % yr; turnover time
RtauFixN                         = [2475, 4125];     % yr; turnover time
% Preferred Values
tauNO3                           = 370;              % yr; turnover time
tauNH4                           = 0.05;             % yr; turnover time
tauPON                           = 0.05;             % yr; turnover time
tauN2                            = 54e3;             % yr; turnover time
tauFixN                          = 33e2;             % yr; turnover time


%% ------------------------------ OXYGEN ----------------------------------
% Bergman et al 2004 COPSE
rp.a.O2.b04.mol                  = 3.7e19;           % mol; atm+ocean oxygen
fp.wthr.c.O2.b04.moloyr          = 0.53e12;          % mol S/yr; modern pyrite oxidative weathering

% from WOCE Atlas online
rr.s.O2.woce.umolar              = [30, 350];        % µmol O2/kg; dissolved surface O2
rr.d.O2.woce.umolar              = [100, 200];       % µmol O2/kg; dissolved deep O2
rr.t.O2.woce.umolar              = rr.s.O2.woce.umolar + rr.d.O2.woce.umolar; 

% Catling et al. 2001 
fp.Hesc.a.H.c01.molyr            = 1e10;             % mol O2/yr; modern day H2 outgassing flux

% Kasting + Canfield, 2012 (Fundamentals of Geobiology, Ch 7)
fp.Hesc.a.H.kc12.molyr           = 9.6e10;           % mol H/yr; modern H escape flux

% from IPCC Ar5 WG1 Ch 6 (Carbon and other Biogeochemical Cycles, 2014)
fp.fire.a.O2.ar14.pgoyr          = 1.61219;          % Pg O2/yr; modern day fire estimated from sec. 6.4.8.1 ~ 1.613 Pg C/yr

% all values without author-quoted ranges are given 25% error bars
rr.a.O2.b04.mol                  = [2.775, 4.625].*1e19; % mol; atm + ocean oxygen 
fr.wthr.c.O2.b04.moloyr          = [0.3975, 0.6625].*1e12;% mol S/yr; modern pyrite oxidative weathering
fr.Hesc.a.H.kc12.molyr           = [7.2 12].*1e10;   % mol H/yr; modern H escape flux
fr.fire.a.O2.ar14.pgoyr          = [1.21, 2.03];     % Pg O2/yr; fire flux estimate


%% ----------------------------- PHOSPHATE --------------------------------
% from Ruttenburg, 2003 (The Global Phosphorus Cycle) 
rr.c.PO4.r03.mol                 = [0.27 1.3].*1e20;   % mol P; crustal rocks, soils and deep marine seds (ie. u + c reservoirs)
rp.c.PO4.r03.mol                 = mean(rr.c.PO4.r03.mol);% mol P; crustal rocks, soils and deep marine seds (ie. u + c reservoirs) - using mean of est range
rp.s.H3PO4.r03.mol               = 87.4.*1e12;         % mol P; total surface ocean dissolved P
rp.d.H3PO4.r03.mol               = 2810.*1e12;         % mol P; total deep ocean dissolved P
rr.t.OP.r03.mol                  = [1.61 4.45].*1e12;  % mol P; total ocean biologic P
fr.wthr.c.PO4.r03.molyr          = [0.245 0.301].*1e12;% mol P/yr; "reactive P" influx (P available for bio-uptake)
fp.wthr.c.PO4.r03.molyr          = mean(fr.wthr.c.PO4.r03.molyr); % mol P/yr; "reactive P" influx (P available for bio-uptake) - mean of est range
fr.burial.z.PO4.r03.molyr        = [0.177 0.242].*1e12;% mol P/yr; burial in ALL marine sediments
fp.burial.n.PO4.r03.molyr        = 1.5e11;             % mol P/yr; burial in continental margin marine seds
fr.burial.n.PO4.r03.molyr        = [0.15 0.223].*1e12; % mol P/yr; burial in continental margin marine seds
fp.burial.z.PO4.r03.molyr        = 1.3e11;             % mol P/yr; burial in abyssal marine seds
fr.burial.z.PO4.r03.molyr        = [0.042 0.13].*1e12; % mol P/yr; burial in abyssal marine seds
fr.burial.z.OP.r03.molyr         = [1.1 4.24].*1e10;   % mol P/yr; burial in ALL marine seds
fr.burial.z.CP.r03.molyr         = [0.4 9.1].*1e10;    % mol P/yr; burial in ALL marine seds
fr.burial.z.FeIIP.r03.molyr      = [0.01 1.43].*1e10;  % mol P/yr; burial in hydrothermal-adjacent marine seds ferric-P
fr.sorb.t.FePO4.r03.molyr        = [0.4 4].*1e10;      % mol P/yr; scavenging/sorption onto Fe-oxides 
fp.burial.z.FePO4.r03.molyr      = 1.5e10;             % mol P/yr; value attributable to burial with CaCO3 ~ FePO4 via sorption sink switching!
fr.ammon.n.H3PO4.r03.molyr       = [0.51 0.84].*1e12;  % mol P/yr; remineralization in coastal seds (modern assumed oxic)
fp.ammon.z.H3PO4.r03.molyr       = 0.41.*1e12;         % mol P/yr; remineralization in abyssal seds
fr.ammon.z.H3PO4.r03.molyr       = [0.287 0.533].*1e12;% mol P/yr; remineralization in abyssal seds (30% uncertainty for preferred est)

% Jahnke 2000 - Earth System Science ch. 14    
rp.c.PO4.j00.mol                 = 0.95.*v.ea.P;     % mol P; crustal abundance of apatite, given 95% crustal P estimated to be hosted in apatite
rr.c.SP.j00.mol                  = [2e-4 1.2e-2].*v.ea.ccrustm;% mol P; apatite in 0.02 - 1.2% of igneous rocks -  assuming most of continental  mass is crystalline, this is an upper estimate
rr.c.CP.j00.mol                  = rp.c.PO4.j00.mol - rr.c.SP.j00.mol;% mol P; apatite in carbonate/sediment rocks,  assuming rest  of apatite is within sediments
rp.c.SP.j00.mol                  = mean(rr.c.SP.j00.mol);% mol P; apatite in igneous rocks  
rp.c.CP.j00.mol                  = mean(rr.c.CP.j00.mol);% mol P; apatite in carbonate/sediment rocks
fp.wthr.c.PO4.j00.mol            = 0.1e12;           % mol P/yr; land to ocean flux (F25,  figure 14-7 and table 14-4) 

% from Paytan + McLaughlin 2007
rp.d.H3PO4.pm07.mol              = 2.9e15;           % mol P; deep ocean P
rp.s.H3PO4.pm07.mol              = 0.1e15;           % mol P; surf ocean P
rp.t.H3PO4.pm07.mol              = 3e15;             % mol P; total ocean
rp.d.H3PO4.pm07.mol              = 2.9e15;           % mol P; deep ocean P
rp.s.H3PO4.pm07.mol              = 0.1e15;           % mol P; surf ocean P
rp.t.H3PO4.pm07.mol              = 3e15;             % mol P; total ocean
rp.c.PO4.pm07.mol                = 0.001.*v.ea.ccrustm./v.const.mm.P; % mol P; crustal P
fr.wthr.c.PO4.pm07.molyr         = 0.75.*[0.6 3].*1e10;% mol P/yr; erosion flux PO4 (modern)
fr.burial.z.PO4.pm07.molyr       = [9.3 34].*0.9e10; % mol P/yr; total reactive P burial in sediments

% Van Capellen and Ingall, 1994
fp.wthr.c.PO4.vci94.molyr        = 3.6e10;           % mol P/yr ; dissolved phosphate erosion
fr.wthr.c.PO4.vci94.molyr        = [0.9 4.5].*1e10;  % mol P/yr ; dissolved phosphate erosion

% from WOCE Atlas online
rr.s.H3PO4.woce.umolar           = [0.05, 2.2];      % µmol P/kg; dissolved phosphorus in surface ocean
rr.d.H3PO4.woce.umolar           = [1.1, 3];         % µmol P/kg; dissolved phosphorus in deep ocean
rr.t.H3PO4.woce.umolar           = rr.s.H3PO4.woce.umolar+rr.d.H3PO4.woce.umolar; % mol P; total ocean phosphate
rp.s.H3PO4.woce.umolar           = mean(rr.s.H3PO4.woce.umolar);
rp.d.H3PO4.woce.umolar           = mean(rr.d.H3PO4.woce.umolar);
rp.t.H3PO4.woce.umolar           = mean(rr.t.H3PO4.woce.umolar);

% all values without author-quoted ranges are given 25% error bars
rr.d.H3PO4.pm07.mol              = [2.175 3.625].*1e15;% mol P; deep ocean P
rr.s.H3PO4.pm07.mol              = [7.5 12.5].*1e13;   % mol P; surf ocean P
rr.t.H3PO4.pm07.mol              = [2.25 3.75].*1e15;  % mol P; total ocean 
rr.c.PO4.pm07.mol                = [3.675 6.125].*1e20;% mol P; crustal P
fr.burial.z.FePO4.r03.molyr      = [1.125 1.875].*1e10;% mol P/yr; value attributable to burial with CaCO3 ~ FePO4 via sorption sink switching!


%% ------------------------------- IRON -----------------------------------
% from Thompson et al. (2019) and references therein
fp.burial.z.FeOH3.t19.tmolyr     = 4.5;              % Tmol Fe/yr; theoretical peak BIF deposition at 2.5 Ga
fp.mantle.t.FeOH3.t19.tmolyr     = 8;                % Tmol Fe/yr; theoretical peak mantle Fe output (8±3 Tmol/yr) at 2.5 Ga
fr.mantle.t.FeOH3.t19.tmolyr     = [5 11];           % Tmol Fe/yr; theoretical peak mantle Fe output in archean
fp.metha.t.CH4.t19.tmolyr        = 3.2;              % Tmol C/yr; theoretical archean methanogenesis 
fp.ferrotrophy.t.CO2.t19.tmolyr  = 31.35;            % Tmol C/yr; mean rate of photo-ferrotrophy primary production, ~1% of total primary production in modern oceans 
fr.ferrotrophy.t.CO2.t19.tmolyr  = [27.6 35.1];      % Tmol C/yr; mean rate of photo-ferrotrophy primary production, ~1% of total primary production in modern oceans 
% rp.d.FeOH3.t19.umolar            = 70;               % µmol Fe/kg; upper constraint on Archean deep ocean [Fe(II)]

% Holland 2006
fp.mantle.t.FeO.h06.molyr        = 3e11;             % mol Fe/yr; modern flux of hydrothermal Fe2+  
% fp.mantle.t.FeO.h06.molyr        = 3e12;             % mol Fe/yr; estimated archean hydrothermal Fe2+ flux (upper limit)

% % Kasting + Canfield, 2012 (Fundamentals of Geobiology, Ch 7)
% fp.burial.z.FeO.kc12.molyr       = 0.9e12;            % mol O2 eq/yr; burial flux of FeO
% fr.burial.z.FeO.kc12.molyr       = [0.5 1.3].*1e12;   % mol O2 eq/yr; burial flux of FeO

% Rudnick + Gao, 2003 (Treatise on Geochemistry)
rp.c.Fe.rg03.mol                 = (0.0504.*v.ea.ccrustm)./v.const.mm.Fe; % mol; upper continental crust iron content ~ 5.04 wt%

% all values without author-quoted ranges are given 25% error bars
fr.burial.z.FeOH3.t19.tmolyr     = [3.375 5.625];    % Tmol Fe/yr; theoretical peak BIF deposition ~ 2.5 Ga
% rp.d.FeOH3.t19.umolar            = [52.5 87.5];      % µmol Fe/kg; Archean deep ocean [Fe(II)]
fr.mantle.t.FeO.h06.molyr        = [2.25 3.75].*1e11;% mol Fe/yr; modern flux of hydrothermal Fe2+  
% fr.mantle.t.FeO.h06.molyr        = [2.25 3.75].*1e12;% mol Fe/yr; estimated archean hydrothermal Fe2+ flux (upper limit)
rr.c.Fe.rg03.mol                 = [0.0409, 0.0726].*v.ea.ccrustm./v.const.mm.Fe; % mol Fe; upper continental crust iron content (4.09 - 7.26 wt%)

%% assign plot symbols to authors
symb.m96   = {'-h','MarkerFaceColor','none'} ;       % == Millero 1996
symb.b95   = {'-p','MarkerFaceColor','none'} ;       % == Bauer 1995
symb.l92   = {'-o','MarkerFaceColor','none'} ;       % == Libes 1992
symb.wod03 = {'-s','MarkerFaceColor','none'} ;       % == Watson + Orr + Doney, 2003
symb.f12   = {'-+','MarkerFaceColor','none'} ;       % == Falkowski 2012
symb.ch18  = {'->','MarkerFaceColor','none'} ;       % == Canuel Hardison 2018
symb.ar14  = {'-^','MarkerFaceColor','none'} ;       % == IPCC Ar 5 2014
symb.wa12  = {'-<','MarkerFaceColor','none'} ;       % == Wallmann Aloisi 2012 (Fundamentals of Geobiology ch 3) 
symb.w08   = {'-d','MarkerFaceColor','none'} ;       % == Wallmann 2008
symb.sm81  = {'-v','MarkerFaceColor','none'} ;       % == Stumm + Morgan 1981
symb.le13  = {'-*','MarkerFaceColor','none'} ;       % == Li + Elderfield 2013
symb.r07   = {'-x','MarkerFaceColor','none'} ;       % == Ridgewell et al 2007
symb.b04   = {'-d','MarkerFaceColor','color'};       % == Bergman et al 2004
symb.w92   = {'-h','MarkerFaceColor','color'};       % == Williams et al. 1992
symb.w12   = {'-p','MarkerFaceColor','color'};       % == Ward 2012
symb.g08   = {'-<','MarkerFaceColor','color'};       % == Gruber 2008
symb.ph18  = {'->','MarkerFaceColor','color'};       % == Palta + Hartnett, 2018
symb.v13   = {'-v','MarkerFaceColor','color'};       % == Voss et al. 2013
symb.pm07  = {'-s','MarkerFaceColor','color'};       % == Paytan + McLaughlin 2007
symb.vci94 = {'-^','MarkerFaceColor','color'};       % == Van Capellen + Ingall 1994
symb.woce  = {'-o','MarkerFaceColor','color'};       % == WOCE online atlas
symb.z12   = {'-h','MarkerFaceColor','c'};           % == Zeebe 2012
symb.f04   = {'-p','MarkerFaceColor','c'};           % == Feely et al. 2004
symb.c01   = {'-s','MarkerFaceColor','c'};           % == Catling et al. 2001
% symb.h78   = {'-d','MarkerFaceColor','c'};           % == Holland 1978 (UNUSED W/O ACTIVE HAZE FLUXES)
symb.h06   = {'-d','MarkerFaceColor','c'};           % == Holland 2006
symb.ms04  = {'-v','MarkerFaceColor','c'};           % == Martin + Sayles, 2004
symb.t19   = {'-^','MarkerFaceColor','c'};           % == Thompson et al. 2019 
symb.kc12  = {'->','MarkerFaceColor','c'};           % == Katling + Canfield, 2012 (Fundamentals of Geobiology ch 7) 
symb.rg03  = {'-<','MarkerFaceColor','c'};           % == Rudnick + Gao, 2003 (Treatise on Geochem)
symb.r03   = {'-o','MarkerFaceColor','c'};           % == Ruttenburg, 2003 (Global Phosphorus cycle)
symb.b03   = {'-V','MarkerFaceColor','y'};           % == Busigny et al. 2003
symb.g09   = {'-o','MarkerFaceColor','y'};           % == Goldblatt et al. 2009
symb.w19   = {'->','MarkerFaceColor','y'};           % == Wong et al. 2019
symb.c17   = {'-<','MarkerFaceColor','y'};           % == Clift 2017
symb.m04   = {'-d','MarkerFaceColor','y'};           % == Mackenzie 2004
symb.j00   = {'-p','MarkerFaceColor','y'};           % == Jahnke 2000
symb.bg14  = {'-^','MarkerFaceColor','y'};           % == Byrne + Goldblatt, 2014
symb.jg15  = {'-s','MarkerFaceColor','y'};           % == Johnson + Goldblatt, 2015

%% Convert units to mol (reservoirs) and mol/yr (fluxes)
% reservoir organization structures house units in 5th field
resstdlist = {rr,rp}; % reservoir ranges and preferred value structures
namelist = {'rr','rp'}; 
for ixr = 1:length(resstdlist)
    rv = resstdlist{ixr}; % either the rp or rr lists
    bxs = fields(rv); % all boxes in that reservoir list
    for ir = 1:length(bxs)
        box = bxs{ir}; % the model box
        sps = fields(rv.(box)); % all species in that box
        for ib = 1:length(sps)
            spe = sps{ib}; % the species
            aus = fields(rv.(box).(spe)); % the author list 
            for is = 1:length(aus)
                auth = aus{is}; % the author
                uz = fields(rv.(box).(spe).(auth)); % all units 
                for ia = 1:length(uz)
                    unit = uz{ia}; % now we get to the juicy bits
                    switch unit
                        case 'umolar' % 1e-6 mol/kg
                            cnv = 1e-6 * v.oc.m.(box); % kg
                        case {'gtc','pgc'} % 1e12 kg C = 1e15 g C == 1 Gt
                            cnv = 1e12 .* (1./v.const.mm.C); % 1e12 kg * mol/kg
                        case {'gtn','pgn'} % 1e12 kg N = 1e15 g N == 1 Gt
                            cnv = 1e12 .* (1./v.const.mm.N); % 1e12 kg * mol/kg
                        case {'gtp','pgp'} % 1e12 kg P = 1e15 g P == 1 Gt
                            cnv = 1e12 .* (1./v.const.mm.P); % 1e12 kg * mol/kg
                        case 'tmol' % 1e12 mol
                            cnv = 1e12; 
                        case 'tgn' % 1e9 kg N = 1e12 g N
                            cnv = 1e9 .* (1./v.const.mm.N); % 1e9 kg * mol/kg
                        case 'tgc' % 1e9 kg C = 1e12 g C
                            cnv = 1e9 .* (1./v.const.mm.C); % 1e9 kg * mol/kg
                        case 'ppb' % ppb = 1e9 mol x/mol atm
                            cnv = 1e-9 .* v.atm.mol; 
                        otherwise % do nothing...
                            cnv = 1; 
                    end
                    % rename into a converted structure!
                    switch namelist{ixr}
                        case 'rr' 
                        crr.(box).(spe).(auth) = {cnv.* rr.(box).(spe).(auth).(unit);symb.(auth)}; 
                        % assign the plot symbol according to the author 
                        case 'rp'
                        crp.(box).(spe).(auth) = cnv.* rp.(box).(spe).(auth).(unit); 
                    end
                end
            end
        end
    end
end
% flux organization structures house units in 6th field 
flxstdlist = {fr,fp}; % flux ranges and preferred value structures
fnamelist = {'fr','fp'}; 

for ixf = 1:length(flxstdlist)
    fv = flxstdlist{ixf}; % either the fp or fr lists
    fxs = fields(fv); % all of the fluxes
    for ir = 1:length(fxs)
        fxx = fxs{ir}; % the flux name
        bxs = fields(fv.(fxx)); % all of the boxes
        for ifx = 1:length(bxs)
            box = bxs{ifx}; % the model box
            sps = fields(fv.(fxx).(box)); % all of the species
            for ib = 1:length(sps)
                spe = sps{ib}; % the species
                aus = fields(fv.(fxx).(box).(spe)); % all of the authors
                for is = 1:length(aus)
                    auth = aus{is}; % the author
                    uns = fields(fv.(fxx).(box).(spe).(auth)); % all of the units
                    for ia = 1:length(uns)
                        unit = uns{ia}; % now we get to the juicy bits
                        switch unit
                            case 'gcm2yr' % grams C/ m^2/yr
                                cnv = 1e-3.*(1./v.const.mm.C).*v.oc.sa; % kg/g * mol/kg * m^2
                            case 'gccm2yr' % grams C/cm^2/yr
                                cnv = 1e-3.*(1./v.const.mm.C).*(v.oc.sa.*100^2); % kg/g * mol/kg * cm^2
                            case 'molm2yr' % mol/m^2/yr
                                cnv = v.oc.sa; 
                            case {'gtcyr','pgcyr'} % 1e12 kg C/yr = 1e15 g C/yr == 1 Gt /yr
                                cnv = 1e12 .* (1./v.const.mm.C); % 1e12 kg * mol/kg
                            case {'gtnyr','pgnyr'} % 1e12 kg N /yr = 1e15 g N/yr == 1 Gt/ yr
                                cnv = 1e12 .* (1./v.const.mm.N); % 1e12 kg * mol/kg
                            case {'gtoyr','pgoyr'} % 1e12 kg O /yr = 1e15 g O/yr == 1 Gt/ yr
                                cnv = 1e12 .* (1./v.const.mm.O); % 1e12 kg * mol/kg
                            case 'tmolyr' % 1e12 mol/yr
                                cnv = 1e12; 
                            case 'tgnyr' % 1e9 kg N/yr = 1e12 g N/yr
                                cnv = 1e9 .* (1./v.const.mm.N); % 1e9 kg * mol/kg
                            case 'tgcyr' % 1e9 kg C/yr = 1e12 g C/yr
                                cnv = 1e9 .* (1./v.const.mm.C); % 1e9 kg * mol/kg
                            otherwise % do nothing...
                                cnv = 1; 
                        end
                        % rename into a converted structure!
                        switch fnamelist{ixf}
                            case 'fr'
                                cfr.(fxx).(box).(spe).(auth) = {cnv.* fr.(fxx).(box).(spe).(auth).(unit);symb.(auth)}; 
                                % assign the plot symbol according to the author
                            case 'fp'
                                cfp.(fxx).(box).(spe).(auth) = cnv.* fp.(fxx).(box).(spe).(auth).(unit); 
                        end

                    end
                end
            end
        end        
    end
end


%% ----------------------------- RESERVOIRS -------------------------------
% preferred reservoir estimates == rx 
rxlst.a.sp = {'O2', 'N2', 'CO2','CH4','NH3','HZ'}; % atm reservoir species
rxlst.a.au = {'b04','w12','f12','ar14','bg14',''}; % atm references
rxlst.s.sp = {'H3PO4','DIC',  'LB',  'fixN', 'N2','NH3','OC','CaCO3','CH4','O2','RN','HNO3','TA','FeOH3'}; 
rxlst.s.au = {'woce', 'f12',  'wa12','woce',  '',  '',   '',  '',     '',   '',  '', 'woce', '',    ''}; 
rxlst.n.sp = {'O2','RN','N2','NH3','HNO3','H3PO4','OC','TA','DIC','CH4','CaCO3','FeOH3'};
rxlst.n.au = {'',   '',  '',  '',   '',    '',    'm04',  '',  '',   '', 'm04',  ''};
rxlst.d.sp = {'H3PO4','DIC', 'fixN', 'RN','N2','TA','HNO3','CaCO3','O2','OC','CH4','FeOH3'};
rxlst.d.au = {'woce', 'f12', 'woce', '',   '',  '',  'woce',  '',   '',  '',  '',    ''}; 
rxlst.z.sp = {'OC', 'O2','RN','N2','HNO3','H3PO4','TA','DIC','CH4','CaCO3','FeOH3'}; 
rxlst.z.au = {'wod03','',   '',  '',  '',    '',     '',  '',  '',  'm04',   ''};
rxlst.u.sp = {'CaCO3','OC','PO4','ON','OP','Fe'}; 
rxlst.u.au = {'w19', 'c17', '','g09',  '',  ''}; 
rxlst.o.sp = {'NH4','CaCO3'}; 
rxlst.o.au = {'g09', 'm04'};
rxlst.c.sp = {'NH4','CaCO3','OC', 'PO4','CP','SP','OP','ON','Fe','N'}; 
rxlst.c.au = {'w12','wa12', 'wa12','r03','j00','j00','f12','f12','rg03','jg15'};
rxlst.m.sp = {'N','C'}; 
rxlst.m.au = {'g09','w19'};
rxlst.t.sp = {'H3PO4','OC','PON','N2', 'HNO3','N',  'DIC', 'CaCO3','RN','TA','O2','CH4','fixN','LB','FeOH3'}; 
rxlst.t.au = {'woce','f12','w12','w12','w12', 'g08','ar14','f12',   '',  '',  '',  '',    'woce', '',   '' }; 
[rx, ~,rxcol] = ReferencePlotSetup(rxlst,crp,'respref',v); 

% range estimates for reservoirs == rxr
rrlst.a.sp = rxlst.a.sp; rrlst.a.au = rxlst.a.au; 
rrlst.s.sp = {'O2',  'RN', 'HNO3','H3PO4','DIC',  'TA', 'LB',   'N2','NH3','OC','CaCO3','CH4','fixN', 'FeOH3'};
rrlst.s.au = {'woce','g08','woce','woce', 'f12',  'm96','wa12',  '',   '',  '',   '',    '',   'woce',  ''}; 
rrlst.d.sp = {'O2',  'RN', 'HNO3','H3PO4','DIC',  'TA','NH3','N2','OC','CH4','CaCO3', 'FeOH3'};
rrlst.d.au = {'woce','g08','woce','woce', 'f12',  'm96','',   '',  '',  '',   '',       ''};
rrlst.n.sp = {'O2','RN','N2','NH3','HNO3','H3PO4','OC','TA','DIC','CH4','CaCO3','FeOH3'};
rrlst.n.au = {'',   '',  '',  '',   '',    '',    'm04', '',  '',   '',  'm04', ''};
rrlst.z.sp = {'DIC','OC','CaCO3','O2','NH3','N2','RN','HNO3','H3PO4','TA','CH4', 'FeOH3'}; 
rrlst.z.au = {'b95', 'm04','m04', '',  '',   '',  '',  '',    '',  	 '', 'wa12', ''}; 
rrlst.u.sp = {'CaCO3','OC', 'PO4','ON','OP','Fe'}; 
rrlst.u.au = {'w19', 'c17','',  'g09', '', 'rg03'}; 
rrlst.o.sp = {'NH4','CaCO3'}; 
rrlst.o.au = {'g09', 'm04'};
rrlst.c.sp = {'NH4','CaCO3','OC', 'PO4','CP','SP','OP','ON','Fe','N'}; 
rrlst.c.au = {'w12','f12', 'f12','r03','j00','j00','f12','f12','rg03','jg15'};
rrlst.m.sp = {'N','C'}; 
rrlst.m.au = {'g09','w19'};
rrlst.t.sp = {'N2', 'HNO3','N',  'DIC', 'PON','H3PO4','OC','RN',  'TA', 'CaCO3','O2',  'CH4','fixN','LB','FeOH3'}; 
rrlst.t.au = {'w12','w12', 'g08','ar14','w12','woce', 'f12','woce','m96','f12',  'woce', '',  'woce','',   ''}; 
[rxr,rsym,rrcol] = ReferencePlotSetup(rrlst,crr,'resrange',v); 

% total carbon/nitrogen in crust/unreactive seds/slab, atm-ocean from sum of all reservoirs
rx.c.C  = rx.c.CaCO3 + rx.c.OC + rx.u.CaCO3 + rx.u.OC + rx.o.CaCO3;
rx.ao.C  = rx.t.CaCO3 + rx.t.CH4 + rx.t.LB + rx.t.DIC + rx.t.OC + rx.a.CO2 + rx.a.CH4; 
rxr.c.C = rxr.c.CaCO3 + rxr.c.OC + rxr.u.CaCO3 + rxr.u.OC + rxr.o.CaCO3;
rxr.ao.C = rxr.t.CaCO3 + rxr.t.CH4 + rxr.t.LB + rxr.t.DIC + rxr.t.OC + rxr.a.CO2 + rxr.a.CH4;
rx.c.N  = rx.c.N + rx.u.ON + rx.o.NH4;
rx.ao.N = rx.t.N + 2.*rx.a.N2 + rx.a.NH3 + rx.t.PON + rx.s.LB./v.const.CNratio;
rxr.c.N = rxr.c.N + rxr.u.ON + rxr.o.NH4;
rxr.ao.N = rxr.t.N + 2.*rxr.a.N2 + rxr.a.NH3 + rxr.t.PON + rxr.s.LB./v.const.CNratio;

% placeholders for values we don't have
bs = {'d','n','z'};
for ib = 1:length(bs)
   rx.(bs{ib}).fixN = rx.(bs{ib}).RN + rx.(bs{ib}).HNO3; 
   rxr.(bs{ib}).fixN = rxr.(bs{ib}).RN + rxr.(bs{ib}).HNO3; 
   if rx.(bs{ib}).fixN == 2e-15
       rx.(bs{ib}).fixN = 1e-15;
   end
   if rxr.(bs{ib}).fixN(1) == 2e-15
       rxr.(bs{ib}).fixN = [1 1].*1e-15;
   end
end

%% -------------------------------- FLUXES --------------------------------
% start with preferred flux values and flux ranges initialized at 0
% for estimates I don't have, using initial values and plotting as white
% symbols, to keep the loop plots working
xsp = {'PO4','NH4','OP','OC','carb','POC','OP','ON','FeOH3'};
for iz = 1:length(xsp)
    fluxx.wthr.(xsp{iz}) = 1e-15; 
    rflux.wthr.(xsp{iz}) = [1e-15 1e-15]; 
    fluxx.meta.(xsp{iz}) = 1e-15;
    rflux.meta.(xsp{iz}) = [1e-15 1e-15];
    fluxx.subduct.(xsp{iz}) = 1e-15;
    rflux.subduct.(xsp{iz}) = [1e-15 1e-15];
    fluxx.volc.(xsp{iz}) = 1e-15;
    rflux.volc.(xsp{iz}) = [1e-15 1e-15];
end 
flus = {'ammon','hdenit','ored','metha'}; sps = {'H3PO4','NH3','DIC','HNO3','OC','RN','N2'}; bxs = {'s','d','n','z'};
for ib = 1:length(bxs)
    bx = bxs{ib};
    for ifu = 1:length(flus)
        fx = flus{ifu}; 
        for is = 1:length(sps)
            sp = sps{is};
            rflux.(fx).(sp).(bx) = [1e-15 1e-15]; 
            fluxx.(fx).(sp).(bx) = 1e-15; 
        end
    end
end
bz = {'s','n','d','z'}; 
for bs = 1:length(bz) 
   rflux.mtrophy.(bz{bs}) = [1e-15 1e-15];
   fluxx.mtrophy.(bz{bs}) = 1e-15;
end
nz = {'n','z'}; sps = {'PO4','CP','ON','OP','OC','POC','CaCO3','FeOH3','FePO4','FeIIP'};
for bs = 1:length(nz)
    for is = 1:length(sps)
        rflux.sed.(sps{is}).(nz{bs}) = [1e-15 1e-15];
        fluxx.sed.(sps{is}).(nz{bs}) = 1e-15;
        rflux.burial.(sps{is}).(nz{bs}) = [1e-15 1e-15];
        fluxx.burial.(sps{is}).(nz{bs}) = 1e-15;
    end
end

sps = {'CaCO3','NH4','SP','CP','OC','ON','OP','FeOH3','Fe2SiO4'}; 
fls = {'subduct','cryst','volc','acc'}; 
for ix = 1:length(sps)
    for il = 1:length(fls)
       fluxx.(fls{il}).(sps{ix}) = 1e-15; 
       rflux.(fls{il}).(sps{ix}) = [1e-15 1e-15]; 
    end
end
dsp = {'RN','NH4','N2','HNO3','NH3','H3PO4','TA','O2','DIC'}; bx = {'z','n'}; 
for id = 1:length(dsp)
    for ib = 1:length(bx)
        rflux.diff.(dsp{id}).(bx{ib}) = [1e-15 1e-15];
    end
end
fluxx.t.oxrmP = 1e-15; 
fluxx.t.axrm = 1e-15; 
fluxx.t.sorb = 1e-15;
fluxx.t.sfw  = 1e-15; 
rflux.t.axrm  = [1e-15 1e-15];
rflux.t.Hesc = [1e-15 1e-15];
rflux.t.sorb = [1e-15 1e-15];
% assign dummy symbols to these references ranges
fxs = fieldnames(rflux);
for is = 1:length(fxs)
    if isa(rflux.(fxs{is}),'double')
         rfsym.(fxs{is}) = {'.','MarkerFacecolor','w'};
    else
        ffs = fieldnames(rflux.(fxs{is})); 
        for iss = 1:length(ffs)
            if isa(rflux.(fxs{is}).(ffs{iss}),'double') == 1
               rfsym.(fxs{is}).(ffs{iss}) = {'.','MarkerFacecolor','w'};
            else
                fffs = fieldnames(rflux.(fxs{is}).(ffs{iss}));
                for isss = 1:length(fffs)
                   if isa(rflux.(fxs{is}).(ffs{iss}).(fffs{isss}),'double') == 1
                       rfsym.(fxs{is}).(ffs{iss}).(fffs{isss}) = {'.','MarkerFacecolor','w'};
                   else
                   end
                end
            end
        end
    end
end

%% now assign the actual flux estimates!
fluxx.gasex.NH3     = cfp.gasex.a.NH3.ar14;  
fluxx.gasex.CO2     = cfp.gasex.a.CO2.wod03;  
fluxx.fix           = cfp.fix.s.N2.ph18;%cfp.fix.s.N2.ar14; % cfp.fix.s.N2.v13; % cfp.fix.s.N2.f12 ;%
fluxx.assim         = cfp.assim.s.NH3.g08; 
fluxx.prod.CO2      = cfp.prod.s.CO2.ch18;   % cfp.prod.s.CO2.wa12;  % cfp.prod.s.CO2.ar14
fluxx.prod.totN     = cfp.fix.s.N2.g08 + cfp.assim.s.NH3.g08; % cfp.prod.s.LB.v13 ;
fluxx.prod.O2       = fluxx.prod.CO2/v.const.COratio; 
fluxx.ferrotrophy.CO2= cfp.ferrotrophy.t.CO2.t19; 
fluxx.ammon.OC.s    = cfp.ammon.s.CO2.ch18; % cfp.ammon.s.CO2.ar14; %
fluxx.ammon.OC.d    = cfp.ammon.d.CO2.ch18; % cfp.ammon.d.CO2.ar14 
fluxx.ammon.OC.n    = cfp.ammon.n.OC.ms04; 
fluxx.ammon.H3PO4.z = cfp.ammon.z.H3PO4.r03; 
fluxx.death         = cfp.death.s.LB.ar14; %flux.death.s.LB.ch18; 
fluxx.export.CaCO3  = cfp.export.s.CaCO3.r07;  
fluxx.export.OC     = cfp.export.s.POC.wod03; % cfp.export.s.POC.r07;  % cfp.export.s.POC.ar14 
fluxx.sed.OC.z      = cfp.sed.z.POC.wod03; %cfp.sed.z.OC.m04;% 
fluxx.sed.CaCO3.z   = cfp.sed.z.CaCO3.m04;
fluxx.burial.PO4.n  = cfp.burial.n.PO4.r03;
fluxx.burial.PO4.z  = cfp.burial.z.PO4.r03;
fluxx.burial.FePO4.z= cfp.burial.z.FePO4.r03;
fluxx.burial.ON.z   = cfp.burial.z.NH4.v13;% cfp.burial.z.NH4.g08; 
fluxx.burial.OC.z   = cfp.burial.z.OC.le13;  % cfp.burial.z.OC.wa12 ;
fluxx.burial.CaCO3.z= cfp.burial.z.CaCO3.wa12; % cfp.burial.z.CaCO3.le13; 
fluxx.burial.CaCO3.n= cfp.burial.n.CaCO3.le13; % cfp.burial.n.CaCO3.ar14; 
fluxx.burial.FeOH3.z= cfp.burial.z.FeOH3.t19; 
fluxx.subduct.NH4   = cfp.subduct.o.NH4.g09; %cfp.subduct.m.NH4.b03; 
fluxx.subduct.ON    = cfp.subduct.u.ON.g09;
fluxx.volc.NH4      = 0.5.*cfp.volc.u.N2.w12; % 0.5.*cfp.volc.u.N2.g08; 
fluxx.volc.CO2      = cfp.volc.u.CO2.ch18;% cfp.volc.u.CO2.ar14;% 
fluxx.wthr.ON       = cfp.wthr.c.ON.v13; % cfp.wthr.c.ON.g08;  
fluxx.wthr.PO4      = cfp.wthr.c.PO4.r03; % cfp.wthr.c.PO4.vci94; 
fluxx.wthr.sil      = cfp.wthr.c.sil.ch18;%cfp.wthr.c.sil.wa12; %cfp.wthr.c.sil.z12;%
fluxx.wthr.carb     = cfp.wthr.c.CaCO3.wa12; % cfp.wthr.c.CaCO3.z12; % cfp.wthr.c.CaCO3.le13; % cfp.wthr.c.CaCO3.wod03; 
fluxx.wthr.oxi      = cfp.wthr.c.OC.b04; % cfp.wthr.c.OC.le13; 
fluxx.wthr.CP       = (crp.c.CP.j00./crp.c.PO4.j00).*fluxx.wthr.PO4; % estimate of total inorganic P weathering given apportionment by Jahnke 2000
fluxx.wthr.SP       = (crp.c.SP.j00./crp.c.PO4.j00).*fluxx.wthr.PO4; % estimate of total inorganic P weathering given apportionment by Jahnke 2000 
fluxx.meta.CaCO3    = cfp.meta.c.CaCO3.wa12; % cfp.meta.c.CaCO3.wod03; %   
fluxx.t.nitr        = cfp.nitr.t.HNO3.g08;  
fluxx.t.denit       = cfp.denit.t.HNO3.g08; % cfp.denit.t.HNO3.b04; % cfp.denit.t.HNO3.ar14; % cfp.denit.t.HNO3.v13; 
fluxx.t.oxrm        = cfp.ammon.s.CO2.ch18+ cfp.ammon.d.CO2.ch18;%cfp.remin.s.CO2.ar14+cfp.remin.d.CO2.ar14; % 
fluxx.t.precip      = cfp.precip.t.CaCO3.wa12; 
fluxx.t.diss        = cfp.diss.t.CaCO3.wa12;  
fluxx.t.revweather  = cfp.revweather.t.TA.le13;%cfp.revweather.t.TA.w08;  
fluxx.t.mantle      = cfp.mantle.t.FeO.h06; 
fluxx.t.methox      = 3.2e13; %IPCC 3rd AR: methane oxidation flux (also == Kasting + Donahue, 1980 via Goldblatt et al. 2008)
fluxx.t.Hesc        = cfp.Hesc.a.H.c01; 
fluxx.t.ncp         = cfp.netcp.t.CO2.f12; %cfp.netcp.t.CO2.ar14;

% flux ranges and symbols!
rflux.gasex.NH3     = cfr.gasex.a.NH3.ar14{1};  
rfsym.gasex.NH3     = cfr.gasex.a.NH3.ar14{2}; 
rflux.gasex.CO2     = cfr.gasex.a.CO2.wod03{1}; % cfr.gasex.a.CO2.ar14; 
rfsym.gasex.CO2     = cfr.gasex.a.CO2.wod03{2}; 
rflux.prod.totN     = cfr.assim.s.NH3.g08{1} + cfr.fix.s.N2.g08{1};  % cfr.prod.s.NH3.v13; 
rfsym.prod.totN     = cfr.assim.s.NH3.g08{2}; 
rflux.prod.CO2      = cfr.prod.s.CO2.ch18{1}; % cfr.prod.s.CO2.wa12;
rfsym.prod.CO2      = cfr.prod.s.CO2.ch18{2}; 
rflux.prod.O2       = rflux.prod.CO2./v.const.COratio;
rfsym.prod.O2       = cfr.prod.s.CO2.ch18{2}; 
rflux.ferrotrophy.CO2= cfr.ferrotrophy.t.CO2.t19{1}; 
rfsym.ferrotrophy.CO2= cfr.ferrotrophy.t.CO2.t19{2}; 
rflux.fix           = cfr.fix.s.N2.ph18{1};%cfr.fix.s.N2.ar14{1};  %cfr.fix.s.N2.v13;  % cfr.fix.s.N2.w12; %cfr.fix.s.N2.g08; 
rfsym.fix           = cfr.fix.s.N2.ph18{2};%cfr.fix.s.N2.ar14{2}; 
rflux.assim         = cfr.assim.s.NH3.g08{1}; 
rfsym.assim         = cfr.assim.s.NH3.g08{2}; 
rflux.ammon.OC.s    = cfr.ammon.s.CO2.ch18{1};  % cfr.remin.s.CO2.ar14{1}; % 
rfsym.ammon.OC.s    = cfr.ammon.s.CO2.ch18{2};
rflux.ammon.RN.s    = cfr.ammon.s.NH3.g08{1}; 
rfsym.ammon.RN.s    = cfr.ammon.s.NH3.g08{2};  
rflux.death         = cfr.death.s.LB.ar14{1}; % cfr.death.s.LB.ch18 
rfsym.death         = cfr.death.s.LB.ar14{2}; 
rflux.export.CaCO3  = cfr.export.s.CaCO3.wa12{1}; % cfr.export.s.CaCO3.r07; 
rfsym.export.CaCO3  = cfr.export.s.CaCO3.wa12{2};
rflux.export.OC     = cfr.export.s.POC.wod03{1}; %cfr.export.s.POC.r07; % cfr.export.s.POC.ar14;
rfsym.export.OC     = cfr.export.s.POC.wod03{2};
rflux.ammon.OC.d    = cfr.ammon.d.CO2.ch18{1}; %cfr.remin.d.CO2.ar14{1}; %
rfsym.ammon.OC.d    = cfr.ammon.d.CO2.ch18{2};
rflux.ammon.RN.d    = cfr.ammon.d.NH3.g08{1}; 
rfsym.ammon.RN.d    = cfr.ammon.d.NH3.g08{2}; 
rflux.ammon.OC.n    = cfr.ammon.n.OC.m04{1}; %cfr.ammon.n.OC.ms04{1};
rfsym.ammon.OC.n    = cfr.ammon.n.OC.m04{2}; %cfr.ammon.n.OC.ms04{2};
rflux.ammon.OC.z    = cfr.ammon.z.OC.m04{1};
rfsym.ammon.OC.z    = cfr.ammon.z.OC.m04{2};
rflux.ammon.H3PO4.z = cfr.ammon.z.H3PO4.r03{1}; 
rfsym.ammon.H3PO4.z = cfr.ammon.z.H3PO4.r03{2};
rflux.diff.DIC.z    = cfr.diff.z.DIC.b95{1};  
rfsym.diff.DIC.z    = cfr.diff.z.DIC.b95{2};  
rflux.sed.OC.z      = cfr.sed.z.OC.wod03{1}; %cfr.sed.z.OC.m04{1}; %   cfr.sed.z.POC.sm81; 
rfsym.sed.OC.z      = cfr.sed.z.OC.wod03{2}; %cfr.sed.z.OC.m04{2}; 
rflux.sed.CaCO3.z   = cfr.sed.z.CaCO3.m04{1};
rfsym.sed.CaCO3.z   = cfr.sed.z.CaCO3.m04{2};
rflux.sed.CaCO3.n   = cfr.sed.n.CaCO3.l92{1}; %cfr.sed.n.CaCO3.m04{1};  
rfsym.sed.CaCO3.n   = cfr.sed.n.CaCO3.l92{2}; %cfr.sed.n.CaCO3.m04{2};  
rflux.burial.ON.z   = cfr.burial.z.NH4.v13{1}; %cfr.burial.z.NH4.g08; 
rfsym.burial.ON.z   = cfr.burial.z.NH4.v13{2};
rflux.burial.OP.z   = cfr.burial.z.OP.r03{1};%cfr.burial.z.PO4.pm07{1};
rfsym.burial.OP.z   = cfr.burial.z.OP.r03{2};%cfr.burial.z.PO4.pm07{2};
rflux.burial.CaCO3.z= cfr.burial.z.CaCO3.m04{1}; %  cfr.burial.z.CaCO3.wa12;% cfr.burial.z.CaCO3.le13; % cfr.burial.z.CaCO3.ar14; 
rfsym.burial.CaCO3.z= cfr.burial.z.CaCO3.m04{2};
rflux.burial.CaCO3.n= cfr.burial.n.CaCO3.le13{1}; %  cfr.burial.n.CaCO3.ar14; 
rfsym.burial.CaCO3.n= cfr.burial.n.CaCO3.le13{2};
rflux.burial.OC.z   = cfr.burial.z.OC.m04{1};%cfr.burial.z.OC.le13{1}; % cfr.burial.z.OC.wa12{1}; 
rfsym.burial.OC.z   = cfr.burial.z.OC.m04{2};%cfr.burial.z.OC.le13{2};
rflux.burial.OC.n   = cfr.burial.n.OC.m04{1};
rfsym.burial.OC.n   = cfr.burial.n.OC.m04{2};
rflux.burial.FeOH3.z= cfr.burial.z.FeOH3.t19{1}; 
rfsym.burial.FeOH3.z= cfr.burial.z.FeOH3.t19{2}; 
rflux.burial.FeIIP.z= cfr.burial.z.FeIIP.r03{1}; 
rfsym.burial.FeIIP.z= cfr.burial.z.FeIIP.r03{2}; 
rflux.burial.FePO4.z= cfr.burial.z.FePO4.r03{1}; 
rfsym.burial.FePO4.z= cfr.burial.z.FePO4.r03{2}; 
rflux.burial.CP.z   = cfr.burial.z.CP.r03{1}; 
rfsym.burial.CP.z   = cfr.burial.z.CP.r03{2}; 
rflux.burial.PO4.z  = cfr.burial.z.PO4.r03{1}; 
rfsym.burial.PO4.z  = cfr.burial.z.PO4.r03{2}; 
rflux.subduct.NH4   = cfr.subduct.o.NH4.g09{1};  %cfr.subduct.m.NH4.b03{1}; 
rfsym.subduct.NH4   = cfr.subduct.o.NH4.g09{2}; %cfr.subduct.m.NH4.b03{2}; 
rflux.subduct.ON    = cfr.subduct.u.ON.g09{1}; 
rfsym.subduct.ON    = cfr.subduct.u.ON.g09{2};  
rflux.subduct.CaCO3 = cfr.subduct.o.CaCO3.w19{1}; 
rfsym.subduct.CaCO3 = cfr.subduct.o.CaCO3.w19{2};
rflux.wthr.ON       = cfr.wthr.c.ON.v13{1}; %  cfr.wthr.c.ON.g08{1}; 
rfsym.wthr.ON       = cfr.wthr.c.ON.v13{2}; 
rflux.wthr.carb     = cfr.wthr.c.CaCO3.wa12{1}; % cfr.wthr.c.CaCO3.z12{1}; % cfr.wthr.c.CaCO3.le13; % cfr.wthr.c.CaCO3.wod03; 
rfsym.wthr.carb     = cfr.wthr.c.CaCO3.wa12{2};
rflux.wthr.sil      = cfr.wthr.c.sil.ch18{1};%cfr.wthr.c.sil.wa12{1}; % cfr.wthr.c.sil.z12{1};
rfsym.wthr.sil      = cfr.wthr.c.sil.ch18{2};% cfr.wthr.c.sil.wa12{2}; %cfr.wthr.c.sil.z12{2};
rflux.wthr.oxi      = cfr.wthr.c.OC.b04{1}; % cfr.wthr.c.OC.le13; 
rfsym.wthr.oxi      = cfr.wthr.c.OC.b04{2};
rflux.wthr.PO4      = cfr.wthr.c.PO4.r03{1}; %cfr.wthr.c.PO4.vci94{1}; %cfr.wthr.c.PO4.pm07;
rfsym.wthr.PO4      = cfr.wthr.c.PO4.r03{2};%cfr.wthr.c.PO4.vci94{2};
rflux.wthr.CP       = (crr.c.CP.j00{1}./crp.c.PO4.j00).*rflux.wthr.PO4; % estimate of total inorganic P weathering given apportionment by Jahnke 2000
rfsym.wthr.CP       = rfsym.wthr.PO4; 
rflux.wthr.SP       = (crr.c.SP.j00{1}./crp.c.PO4.j00).*rflux.wthr.PO4; % estimate of total inorganic P weathering given apportionment by Jahnke 2000 
rfsym.wthr.SP       = rfsym.wthr.PO4; 
rflux.volc.NH4      = 0.5.*cfr.volc.u.N2.w12{1}; % 0.5.*cfr.volc.u.N2.g08;  
rfsym.volc.NH4      = cfr.volc.u.N2.w12{2};
rflux.volc.CO2      = cfr.volc.u.CO2.w19{1}; %cfr.volc.u.CO2.ch18{1}; %cfr.volc.u.CO2.wa12{1}; % cfr.volc.u.CO2.ar14; % 
rfsym.volc.CO2      = cfr.volc.u.CO2.w19{2};%cfr.volc.u.CO2.ch18{2}; %cfr.volc.u.CO2.wa12{2};
rflux.volc.carb     = cfr.volc.u.carb.w19{1};
rfsym.volc.carb     = cfr.volc.u.carb.w19{2};
rflux.volc.OC       = cfr.volc.u.OC.w19{1}; %cfr.volc.u.ON.g09{1}.*v.const.CNratio; 
rfsym.volc.OC       = cfr.volc.u.OC.w19{2}; %cfr.volc.u.ON.g09{2};
rflux.meta.CaCO3    = cfr.meta.c.CaCO3.wa12{1}; % cfr.meta.c.CO2.wod03; %
rfsym.meta.CaCO3    = cfr.meta.c.CaCO3.wa12{2};
rflux.meta.OC       = cfr.meta.c.OC.wa12{1}; % cfr.meta.c.OC.b04; 
rfsym.meta.OC       = cfr.meta.c.OC.wa12{2};
rflux.t.oxrm        = rflux.ammon.OC.d+rflux.ammon.OC.s;
rfsym.t.oxrm        = rfsym.ammon.OC.d; 
rflux.t.precip      = cfr.precip.t.CaCO3.wa12{1}; %cfr.precip.t.CaCO3.m04{1}; %
rfsym.t.precip      = cfr.precip.t.CaCO3.wa12{2}; %cfr.precip.t.CaCO3.m04{2};%
rflux.t.diss        = cfr.diss.t.CaCO3.wa12{1};%cfr.diss.t.CaCO3.f04{1}; % 
rfsym.t.diss        = cfr.diss.t.CaCO3.wa12{2}; %cfr.diss.t.CaCO3.f04{2};%
rflux.t.nitr        = cfr.nitr.t.HNO3.g08{1}; %cfr.nitr.t.HNO3.b04;
rfsym.t.nitr        = cfr.nitr.t.HNO3.g08{2};
rflux.t.denit       = cfr.denit.t.HNO3.g08{1}; %cfr.denit.t.HNO3.b04; %cfr.denit.t.HNO3.v13; 
rfsym.t.denit       = cfr.denit.t.HNO3.g08{2};
rflux.t.revweather  = cfr.revweather.t.TA.le13{1};%cfr.revweather.t.TA.w08{1}; 
rfsym.t.revweather  = cfr.revweather.t.TA.le13{2};%cfr.revweather.t.TA.w08{2}; 
rflux.t.mantle      = cfr.mantle.t.FeO.h06{1}; %cfr.mantle.t.CH4.b04{1}; 
rfsym.t.mantle      = cfr.mantle.t.FeO.h06{2}; %cfr.mantle.t.CH4.b04{2}; 
rflux.t.methox      = cfr.methox.t.CH4.ar14{1};
rfsym.t.methox      = cfr.methox.t.CH4.ar14{2};
rflux.t.sorb        = cfr.sorb.t.FePO4.r03{1};
rfsym.t.sorb        = cfr.sorb.t.FePO4.r03{2};
rflux.t.ncp         = cfr.netcp.t.CO2.f12{1}; %cfr.netcp.t.CO2.ar14{1}; % cfr.netcp.t.CO2.ch18{1}; 
rfsym.t.ncp         = cfr.netcp.t.CO2.f12{2};
rflux.t.sfw         = cfr.sfw.t.CaCO3.w08{1};
rfsym.t.sfw         = cfr.sfw.t.CaCO3.w08{2};

%% what are some implied timescales from these literature estimates?
fs = {0.5.*rflux.fix, rflux.assim, rflux.death, rflux.t.nitr, rflux.t.denit};
rs = {rxr.s.fixN, rxr.s.fixN, rxr.s.LB, rxr.t.HNO3, rxr.t.HNO3};

for ifs = 1:length(fs)
   [itau{ifs}] = ImpliedTimescale(fs{ifs},rs{ifs}); 
end

% and get some different species estimates for the same kind of flux
ffs = {'meta','volc','wthr'}; 
sfs{1} = {'CaCO3','OC'};       rfs{1} = {rxr.c.CaCO3,rxr.c.OC}; 
sfs{2} = {'carb','OC'};        rfs{2} = {rxr.u.CaCO3, rxr.u.OC}; 
sfs{3} = {'oxi','carb','PO4'}; rfs{3} = {rxr.c.OC, rxr.c.CaCO3, rxr.c.PO4};

for irs = 1:length(ffs)
    fls = ffs{irs}; spes = sfs{irs}; 
    for is = 1:length(spes) 
        [iitau.(fls).(spes{is})] = ImpliedTimescale(rflux.(fls).(spes{is}),rfs{irs}{is}); 
    end
end

%% subfunction for setting up plot details with structure inputs!
function [values,symbol,color] = ReferencePlotSetup(list,ref,flag,v)
% use a structure input that includes all reservoirs/fluxes, species, and
% author labels to set up the plot reference values, color, and symbol
switch flag
    % recall that the structure of resrvoirs is :
    % (LISTNAME).RESERVOIR.SPECIES.AUTHOR = VALUE
    case {'resrange', 'respref'} % reservoirs!
        if strcmp(flag,'resrange')
            norefval = [1e-15 1e-15]; % use a range of out-of-bounds values
        else
            norefval = 1e-15; % use a single out-of-bounds value
        end
        xss = fieldnames(list); % the values structure
        for irs = 1:length(xss)
            spx = list.(xss{irs}).sp; % reservoir's species
            for ixs = 1:length(spx)
                spe = spx{ixs}; % the species
                if isfield(ref,xss{irs})  % make sure the reservoir is in the reference list!
                    aux = list.(xss{irs}).au; % author references for that reservoir
                    % check if the reservoir has a reference for this species
                    if isfield(ref.(xss{irs}),spe) 
                        if length(ref.(xss{irs}).(spe).(aux{ixs})) == 2
                            % assign value to the reservoir output structure and plot symbol based on author
                            values.(xss{irs}).(spe) = ref.(xss{irs}).(spe).(aux{ixs}){1}; 
                            symbol.(xss{irs}).(spe) = ref.(xss{irs}).(spe).(aux{ixs}){2};
                        else
                            values.(xss{irs}).(spe) = ref.(xss{irs}).(spe).(aux{ixs}); 
                            symbol.(xss{irs}).(spe) = {'.','MarkerFaceColor','none'};
                        end
                    else  % for those reservoirs that don't have references, assign a dummy symbol and value!
                        symbol.(xss{irs}).(spe) = {'.','MarkerFaceColor','none'};
                        values.(xss{irs}).(spe) = norefval;
                    end
                else  % for those reservoirs that don't have references, assign a dummy symbol and value!
                    symbol.(xss{irs}).(spe) = {'.','MarkerFaceColor','none'};
                    values.(xss{irs}).(spe) = norefval;
                end
                [color.(xss{irs}).(spe)] = ReferenceColorSetup(spe,v); % assign a plot color according to the species
                % if the plot symbol includes the 'color' designation,
                % replace this placeholder with the actual color!
                symcolor = symbol.(xss{irs}).(spe){3}; 
                switch symcolor
                    case 'color'
                        symbol.(xss{irs}).(spe){3} = color.(xss{irs}).(spe); 
                    otherwise % do nothing! the color designation otherwise is 'none'
                end
                        
            end
        end
    % recall that the structure of fluxes is :
    % (LISTNAME).FLUX.RESERVOIR.SPECIES.AUTHOR = {VALUE,SYMBOL}
    case 'fluxes'
        xss = fieldnames(list); % the reservoirs in the flux structure
        for irs = 1:length(xss)
            rsx = xss{irs}; % the reservoir name
            fxx = list.(rsx).flx; % the fluxes in each reservoir
            spx = list.(rsx).sp; % fluxed species in each reservoir
            aux = list.(rsx).au; % the author list in each reservoir
            for ixs = 1:length(spx)  % ixs = the length of all the lists!
                spe = spx{ixs}; % the species
                flu = fxx{ixs}; % the flux
                aus = aux{ixs}; % the author we want to use
                if isfield(ref,flu) % make sure that the flux is actually included in the list
                    if isfield(ref.(flu).(rsx),spe) % make sure if the reference list has a flux for that species in that reservoir
                        % check if the flux has a reference for this species
                        if isfield(ref.(flu).(rsx).(spe),aux{ixs}) 
                            if length(ref.(flu).(rsx).(spe).(aux{ixs})) == 2
                                % assign value {1} to the reservoir output structure and plot symbol {2} based on author
                                values.(flu).(rsx).(spe) = ref.(flu).(rsx).(spe).(aus){1}; 
                                symbol.(flu).(rsx).(spe) = ref.(flu).(rsx).(spe).(aus){2};
                            else
                                values.(flu).(rsx).(spe) = ref.(flu).(rsx).(spe).(aus); 
                                symbol.(flu).(rsx).(spe) = {'.','MarkerFaceColor','none'};
                            end
                        else  % for those fluxes that don't have references, assign a dummy symbol and value!
                            values.(flu).(rsx).(spe) = [1e-15 1e-15];
                            symbol.(flu).(rsx).(spe) = {'.','MarkerFaceColor','none'};
                        end 
                    else  % for those fluxes that don't an entry in the reference list, assign a dummy symbol and value!
                        values.(flu).(rsx).(spe) = [1e-15 1e-15];
                        symbol.(flu).(rsx).(spe) = {'.','MarkerFaceColor','none'};
                    end
                else % make a dummy symbol and value
                    values.(flu).(rsx).(spe) = [1e-15 1e-15];
                    symbol.(flu).(rsx).(spe) = {'.','MarkerFaceColor','none'};
                end 

               [color.(flu).(rsx).(spe)] = ReferenceColorSetup(spe,v); 
                % if the plot symbol includes the 'color' designation,
                % replace this placeholder with the actual color!
                symcolor = symbol.(flu).(rsx).(spe){3}; 
                switch symcolor
                    case 'color'
                        symbol.(rsx).(spe){3} = color.(flu).(rsx).(spe); 
                    otherwise % do nothing! the color designation otherwise is 'none'
                end
            end
        end
end
end

%% subfunction for assigning colors to references, based on elements
function [colr] = ReferenceColorSetup(spe,v)
es = {'H','C','O','N','P','Fe'}; % the elements we care about
switch spe
    case {'OC','POC','ON','OP','DIC','RN','fixN','PON','CP','SP'} % in this case, ignore the next eval and just assign the final element as the key color variable
        eidx = length(spe); 
    case 'H3PO4'
        eidx = 3; % use P!
    case 'HNO3'
        eidx = 2; 
    case {'FeOH3','Fe'}
        eidx = 1:2; 
    otherwise
        eidx = 1; 
end
for iex = 1:length(es) % iterate through the element list
    if contains(spe(eidx),es{iex}) == 1 % look at the element in the index provided
        colr = v.color.(es{iex}) ; % the species color = the key element color
    elseif strcmp(spe,'TA') || strcmp(spe,'LB') || strcmp(spe,'DB') || strcmp(spe,'HZ') || strcmp(spe,'sil') || strcmp(spe,'org')% use C as color
        colr = v.color.C;
    else % try the next element in the list
    end
end 

end


%% calculate implied timescales
function [impliedtau] = ImpliedTimescale(flux,reservoir)
% implied turnover timescales are calculated by taking the reservoir of the
% product and dividing by the flux into or out of the reservoir (ie. burial
% timescale for N == buried N reservoir / est. flux of burial)

impliedtau = reservoir ./ flux; 

end


