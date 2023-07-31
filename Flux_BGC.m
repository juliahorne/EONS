%% ============================= EONS Model =============================== 
% Julia Horne, 2018
%
% Equations for all of the biogeochemical fluxes in the model, including: 
% NH3 photolysis, continental erosion, primary production, death of 
% living biomass, export production, remineralization, burial, CaCO3 
% precipitation and dissolution, nitrification/denitrification, etc.

% Comments regarding sources and stoichiometry accompany all model fluxes- 
% some more detailed comments are included as footnotes at the end of this
% function (denoted by *** footnote X ***)

function flux = Flux_BGC(mt,conc,r,T,inp,v)

% preliminary details
refsize = size(r.a.N2);                                     % reference size of all flux/reservoir vectors! used in preallocation
[vr]    = SetStoichiometricRatios('Canon');                 % what organic matter stoichiometry we will use
[t]     = CheckTheTime(mt,inp);                             % make sure we use real time in a full model run
fO2     = r.a.O2./v.atm.O2pal;                              % PAL; relative abundance of oxygen in atmosphere
tdep    = Forcings(t,v);                                    % time-dependent forcings (tdep)
[tau]   = Timescales(tdep,inp,v);                           % flux rate controls 
c       = VolumetricConcentrations(r,conc,v);               % calculate volumetric concentrations for all species
lim     = Limitations(r,conc,c,v);                          % functional limitations on biosphere

%%                            Boundary forcings
% Erosion has non-linear relationship with topography (Montgomery + Brandon
% 2002) and should have near-modern rates by mesoarchean; we don't include 
% topography, so instead related to continental emergence and ensure 
% sedimentation is within 2 mag of modern day - this gets a curve that's
% around 30% modern erosion by 3 Ga
xcont = v.ea.Si ./ r.c.SiO3;                                % 1/continental fraction (ie. starts at ~20, evolves to ~1)
xeros = xcont .^ 1.5;                                       % relationship to erosion (continental emergence) 
% OPTIONAL: modify accretion rate (for both oceanic crust and sediments) and/or productive shelf area according to continental size
f_acc.o = v.f.acc;% ./ xcont;                                 % accretion fraction w/r/to continental size (increasing)
f_Pshelf = v.f.Pshelf;% ./ xcont;                             % productive shelf area w/r/to continental size (increasing)

% ------------------------ Mantle Reductant Influx ------------------------
% A key redox boundary condition: FeO comes from the mantle and goes into 
% the deep ocean, where it is mixed and either oxidized by free O2 or via 
% anaerobic photosynthesis. Baseline assumption is that this decreases from
% imp.fman x modern to modern flux (3e11 mol/yr; Holland, 2002) as mantle 
% outgassing eases (Sleep + Zahnle, 2001)
flux.mantle.FeO = EvolveParam(t,4e9,inp.fman.*3e11,3e11,refsize); % mol FeO/yr; influx linearly decreasing through time 

if isfield(inp,'testname')                                  % if running a mantle sensitivity test, this flux may be constant
    switch inp.testname
        case 'constant'
            flux.mantle.FeO = inp.fman .* 3e11 .* ones(refsize); % mol FeO/yr; influx constant
    end
end
% --------------------------- Mantle Outgassing ---------------------------
% We use evolving mantle reservoirs to determine the outfluxes of carbon,
% nitrogen, silicate-bound P, SiO3 and Fe to the surface, assuming only 10%
% of the upper mantle will "degas" at any given time (f_UMgas = 0.1). 
% Silicates are emplaced as LIPs (ie. SiO3, SP, and Fe2SiO4)

els = {'C','N','P','SiO3','Fe'}; 
for ie = 1:length(els)
    msp = els{ie}; mn = msp;  
    switch msp
        case 'Fe' % treat this as a flux in units of mol Fe2SiO4/yr, not mol Fe/yr
            mn = 'Fe2SiO4'; xmol = 0.5;
        otherwise 
            xmol = 1; 
    end
    flux.mantle.(mn) = xmol .* v.f.UMgas .* (r.m.(msp) ./ tau.mantle);% mol/yr; outflux from mantle
end
% Carbon and nitrogen are outgassed as volatiles (CO2/CH4, N2/NH3), in 100:1 
% ratios oxidized/neutral to reduced; estimated from 2 orders of magnitude
% difference in modern CO2:CH4 volcanic outgassing ratios (Wong et al.
% 2019)
flux.mantle.CO2 = (1 - v.f.UMred) .* flux.mantle.C;         % mol CO2/yr; 
flux.mantle.CH4 = v.f.UMred .* flux.mantle.C;               % mol CH4/yr;
flux.mantle.N2  = 0.5 .* (1 - v.f.UMred) .* flux.mantle.N;  % mol N2/yr;
flux.mantle.NH3 = v.f.UMred .* flux.mantle.N;               % mol NH3/yr; 

%%                         Atmosphere reactions
%------------------- Ammonia + Methane Photo-oxidation --------------------
% Ammonia is soluble in the ocean as well as being irreversibly
% photolysized in the atmosphere, where high energy photons break up the
% molecule, and react hydrogen and available oxygen:
%        NH3 --3hv--> 1/2 N2 + 3H+ + 3/4 O2 -> 1/2 N2 + 3/2 H2O 
% This occurs at a rate directly dependent on oxygen level (faster rxn with 
% more O2), and the reaction itself depends on [NH3]. We dub this 'ammox' 
% for ammonia oxidation, which takes the form: 
%                          Flux = kamx[O2] rNH3
flux.ammox = zeros(refsize);                                % initialize vector
kamx = 1 ./ (tau.ammox .* 0.21);                            % /yr [O2]; oxygen concentration sensitive rate constant, with respect to modern oxygen concentration
O = r.a.O2./v.atm.mol;                                      % abundance O2; atmospheric oxygen concentration
for io = 1:length(t)
    if (r.a.O2(io) > 0) && (r.a.NH3(io) > 0)                % ensure this goes to zero if NH3 or O2 are depleted and not negative
        flux.ammox(io) = kamx .* O(io) .* r.a.NH3(io);      % mol N/yr; NH3 photo-oxidation flux  
    else
        flux.ammox(io) = 0;
    end
end
% Ammonia is broken down by UV rays even when oxygen is not available.
% In the anoxic atmosphere the reaction involves CO2 and CH4:
%              8NH3 + 4CO2 --hv--> 4N2 + 4CH4 + 4H2O + 2O2 
%                    CH4 + 2O2 --hv--> CO2 + 2H2O
% With a net reaction:
%                  8NH3 + 3CO2 --hv--> 4N2 + 3CH4 + 6H2O
% We dub this 'photolysis' (pholys). This flux occurs over a longer
% timescale for higher ammonia concentrations, inndependent of [O2]: 
%                           Flux = kphy rNH3
flux.pholys = zeros(refsize);                               % initialize vectors
kphy = 1 ./ tau.pholys;                                     % /yr; longer anoxic NH3 lifetime (Kasting, 1982)
for io = 1:length(t)
    if r.a.NH3(io) > 0                                      % ensure this goes to 0 if NH3 is depleted and not negative
        flux.pholys(io) = kphy .* r.a.NH3(io);              % mol N/yr; NH3 photolysis flux  
    else
        flux.pholys(io) = 0;
    end
end
% Oxidation of CH4 to CO2 requires free oxygen in the atmosphere and
% occurs at a rate proportional to atm OH- levels (which are not considered 
% here; Goldblatt et al. 2006)
%                       2O2 + CH4 --4hv--> CO2 + 2H2O 
% Oxidation rate is dependant on an effective rate constant (Keff) which
% has been approximated by Colin based on ATMOS model output from Daniel 
% Garduno-Ruiz (2023). Generally:
%                     Flux = Keff([CH4][O2])^1/2 * L_O2 
% We impose a limitation on this flux at O2 mixing ratios below 1e-12;
% otherwise this parameterization is inclined to crash the oxygen reservoir
% when it turns on.
M = r.a.CH4./v.atm.mol;                                     % abundance CH4; atmospheric methane concentration
Keff = zeros(size(M)); flux.methox = zeros(size(M));  Ol = zeros(size(M)); 
for ikf = 1:length(r.a.O2) 
    Keff(ikf) = 3e20.*4e3.^((-O(ikf).^0.4)./(O(ikf) + 1e-4).^0.4);
    Ol(ikf)   = O(ikf) ./  (O(ikf) + 1e-12);                % smoothly increase the flux as oxygen rises, so it doesn't crash the reservoir
    if (r.a.CH4(ikf) > 0) && (r.a.O2(ikf) > 0)
        flux.methox(ikf) = Keff(ikf).*(M(ikf).*O(ikf)).^0.5 .* Ol(ikf);% mol O2/yr; oxidation flux 
    else
        flux.methox(ikf) = 0; 
    end
end

%---------------------------- Hydrogen Escape -----------------------------
% Hydrogen loss at the upper atmosphere acts as a net oxidant, and supply
% is limited by diffusion from the lower atmosphere. Methane and ammonia
% photolysis supply H, and thus total hydrogen loss is proportional to the
% total atmospheric reservoir of those gases, and to a rate controlled by 
% the molar atmospheric volume, assumed constant (Goldblatt et al., 2006)
%                       NH3 + 3hv -> N3- + 3H+ (out)
%                       CH4 + 4hv -> C4- + 4H+ (out) 
%                       C4- + H2O -> 1/2CH4 + 1/2CO2
% The flux of hydrogen escaping is determined from the mixing ratios (f_H) 
% of each gas with H and a constant escape flux estimated in Walker (1977).
% Formulation is from Claire et al. (2006):
%                    k_esc = 2.5e13 (H2 molecules/cm2 s)
%                   Flux = 2k_esc * f_H (H molecules/yr)
kesc = 2*v.atm.kesc*v.conv.atmunits;                        % mol H/yr; hydrogen escape flux per mol H
flux.Hesc.CH4 = zeros(refsize);flux.Hesc.NH3 = zeros(refsize); fg = zeros(refsize);
gs = {'CH4','NH3'}; frac = {(1/4), (1/3)}; 
for ig = 1:length(gs)
    for ir = 1:length(r.a.(gs{ig}))
        fg(ir) = r.a.(gs{ig})(ir)./v.atm.mol;               % CH4, NH3 mixing ratio
        if r.a.(gs{ig})(ir) > 0 
            flux.Hesc.(gs{ig})(ir) = fg(ir) .* kesc .* frac{ig};% mol gas/yr; gas loss from H escape
        else
            flux.Hesc.(gs{ig})(ir) = 0; 
        end
    end
end
flux.Hesc.H = 3.*flux.Hesc.NH3 + 4.*flux.Hesc.CH4;          % mol H/yr; total H escape

%%                         Subaerial Weathering
%-------------- Carbonate and Silicate Weathering Modifiers ---------------
% Weathering fluxes are sensitive to atmospheric pCO2 and temperature, as 
% combination with rainwater generates weak carbonic acid and eventually 
% forms soluble bicarbonate ions, transporting alkalinity/DIC to the ocean 
% in a climate-stabilizing feedback (Walker et al. 1981; Brady + Carrol, 
% 1994; Bergman et al. 2004)
[ws] = WeatheringSensitivities(t,T,r.a.CO2,conc.z.CO2,tdep,v);

%------------------------------ Weathering --------------------------------
% Oxidative weathering of organic carbon, operates essentially as 
% remineralization on the continent by atmospheric O2 (Bergman et al 2004):
%                    CH2O + O2(a) --> CO2 + H2O
% This releases NH3 and H3PO4 in Redfield ratio (RR). We calculate out of 
% this ratio because P can be scavenged in anoxic seds and may not be 
% buried strictly along RR!
flux.wthr.oxi = zeros(refsize); flux.wthr.ON = zeros(refsize); flux.wthr.OP = zeros(refsize); 
for io = 1:length(r.a.O2)
    sps = {'OC','OP','ON'};
    for is = 1:length(sps)
        if r.a.O2(io) > 0 && r.c.(sps{is})(io) > 0
            if strcmp(sps{is},'OC')
            flux.wthr.oxi(io) = sqrt(fO2(io))  .* (r.c.(sps{is})(io)...
                ./ tau.woxi);                               % mol C,O2/yr
            else
            flux.wthr.(sps{is})(io) = sqrt(fO2(io)) .* (r.c.(sps{is})(io)...
                ./ tau.woxi);                               % mol N,P/yr

            end
        end
    end
end

% Inorganic carbon (CaCO3) weathering on continents uses atmospheric CO2 
% and water, producing alkalinity:
%             CaCO3(c) + CO2(a) + H2O --> Ca^2+ + 2HCO3^- 
flux.wthr.carb = ws.carb .* (r.c.CaCO3 ./ tau.wcarb);       % mol C/yr
% Silicate weathering uses CO2 and water from the atmosphere to form weak 
% carbonic acid, forming clays and bicarbonate:
%             CaSiO3 + 2CO2(a) + H2O --> 2HCO3^- + Ca^2+ + SiO2
flux.wthr.sil = 2 .* ws.sil .* (r.c.SiO3 ./ tau.wsil);      % mol C/yr; atm CO2 used per mol of silicate (CaSiO3)

% Carbonate-bound phosphate (apatite = CP) is removed via carbonate 
% weathering:
%      Ca5(PO4)3(OH) + CO2(a) + H2O -> 5Ca^2+ + 3PO4 + HCO3^- + H2O
% we simplify to:
%       Ca3(PO4)2 + 3CO2(a) + 3H2O(a) -> 3Ca^2+ + 2H3PO4 + 3CO3^2-
flux.wthr.CP = ws.carb .* (r.c.CP ./ tau.wcarbP);           % mol P/yr
% The "primordial" phosphorus source on continents is a silicate-bound P 
% reservoir sourced from mantle LIP eruption, replenished by 
% recrystallization of OP + CP sediments and the burial of Fe-oxide bound P
% (FePO4). Silicate-bound P is weathered with CO2 to produce first apatite 
% and then dissolved phosphate and alkalinity:
%      P2O5 + 3CaSiO3 + 6CO2 + 6H2O -> Ca3(PO4)2 + 3SiO2 + 6H2CO3
%                   -> 3Ca^2+ + 3SiO2 + 2H3PO4 + 6HCO3^-
% Simplified to reduce untracked species:
%             PO4(c) + 3CO2(a) + 3H2O(a) -> H3PO4 + 3HCO3^- 
flux.wthr.SP = ws.sil .* (r.c.SP ./ tau.wsilP);             % mol P/yr
% Total inorganic P weathering is assumed to be ~ 10e10 mol/yr, split 
% between organic and carbonate + silicate (inorganic) sources:
flux.wthr.PO4 = flux.wthr.CP + flux.wthr.SP; 
% Inorganic P weathering (CP and SP) increases when fungi evolve and begin
% selective remineralization of this species on the continent. Therefore,
% the timescale of inorganic P weathering is a forcing in the model.

% Silicate-bound ammonium is removed from buddingtonite via silicate 
% weathering:
% 4(NH4AlSi3O8) + 4CO2(a) + 6H2O -> 4NH4^+ + 2Al2Si2O5(OH)4 + 8SiO2 + 4HCO3^-
% we simplify this to a stoichiometry similar to regular silicate 
% weathering to reduce untracked species:
%          (NH4)2SiO3 + 2CO2(a) + H2O(a) --> 2NH4^+ + SiO2 + 2HCO3^- 
flux.wthr.NH4 = ws.sil .* (r.c.NH4 ./ tau.wsil);            % mol N/yr

% Reduced iron Fe(II) is a large component of initial continental
% materials (fayalite olivine), and the weathering of this reduced 
% iron draws CO2 and O2 to create oxidized hematite, via:
%      Fe2SiO4 + CO2 + H2O + (1/2)O2 --> Fe2O3 + SiO2 + H2CO3
% Oxidized iron remains on the continent, while produced DIC is transported
% to the ocean.
flux.wthr.Feox = zeros(refsize);  
for it = 1:length(t)
    if r.c.Fe2SiO4(it) > 0 && r.a.O2(it) > 0 
        flux.wthr.Feox(it) = sqrt(fO2(it)) .* ws.sil(it)...
            .* (r.c.Fe2SiO4(it) ./ tau.wsil);               % mol Fe/yr
    end
end

%%                          Primary Production
%---------------------- Fixation and Assimilation -------------------------
% initialize vectors  
flux.assim = zeros(refsize); flux.ferrotrophy.FeO = zeros(refsize); 
flux.fixation.FeO = zeros(refsize); flux.fix.newN = zeros(refsize);
fRN = r.s.RN./(r.s.RN + r.s.HNO3);                          % ratio of reduced N to total bio N

for iT = 1:length(t)

% --------------------------- Photosynthesis ------------------------------
% Producers generate organic matter from C, N, and P nutrients in the
% presence of light; this reaction is primarily limited by phosphate. The
% ratio of reactants to products depends on the stoichiometry of OM and 
% fixed N.
% Stoichiometry from Gruber (2008):
%    106CO2 + 16(NO3^-, NH4^+) + HPO4^2- + (78,48)H2O + (18H^+, 14OH^-)
%                --> C_106 H_175 O_42 N_16 P + (150,118)O2
% Canonical Stoichiometry (Fundamentals of Geobiology, Konhauser et al. + Berman-Frank et al. 2008):
%          106CO2 + 16(NH3, NH4^+, HNO3) + H3PO4 + (106,106,122)H2O 
%         --> C_106 H_263 O_110 N_16 P + (106,106,138)O2 + (0,16,0)H^+
% NOTE: even if charged species NH4+ is taken as the fixed N source, the 
% reaction is alkalinity neutral because of the maintained OM composition
    if r.s.H3PO4(iT) > 0
        flux.assim(iT)= r.s.H3PO4(iT) .* lim.C(iT) .* lim.N(iT) .* lim.RS(iT)...
            .* v.const.CPratio .* (1/tau.assim);            % mol C/yr
    end
    % oxygenic photosynthesis evolves after 1 billion years
    flux.assim(iT) = flux.assim(iT) .* tdep.photo(iT);
    
% --------------------------- Photoferrotrophy ----------------------------    
% This production pathway is limited by light, ferrous iron availability,
% and limited to anoxic settings (Crowe et al. 2008). Due to dependence on 
% anoxia, we assume that it's fixed N can only be either ammonia/ammonium, 
% not nitrate. Stoichiometry modified from Goldblatt et al. (2006) to be 
% alkalinity neutral: 
%            4FeO + CO2 + 11H2O --hv--> 4Fe(OH)3 + CH2O + 4H2O
% ie.      4FeO + CO2 + 16/106(NH3, NH4^+) + 1/106 H3PO4 + 11H2O
%   --> 4Fe(OH)3 + 1/106 (CH2O)_106 (NH3)_16 (H3PO4) + 4H2O + (0,16/106)H^+
    if r.s.FeO(iT) > 0
        flux.ferrotrophy.FeO(iT) = r.s.FeO(iT) .* lim.C(iT) .* lim.P(iT) .* lim.RN(iT)...
            .* lim.RS(iT) .* (1 - lim.Phx(iT)) .* (1/tau.assim);% mol Fe/yr
    end

%------------------------------ N2 Fixing  --------------------------------
% Newly fixed nitrogen has to be generated when there is plentiful P but 
% insufficient fixed nitrogen; this is energy intensive due to N2's triple 
% bond, so power is sourced by transferring ATP to ADP. *** footnote 1 ***
% Modern (oxygenic) fixation by cyanobacteria:
%                      N2 + 3H2O --> 2NH3 + 3/2 O2
% This is essentially a photosynthetic reaction, leading to immediate
% uptake of newly fixed nitrogen; the entire reaction is therefore:
%   106CO2 + 8N2 + H3PO4 + 130H2O --> (CH2O)_106(NH3)_16(H3PO4) + 118O2
% 
% oxygenic nitrogen fixation == 'fix'
    if r.s.H3PO4(iT) > 0
        flux.fix.newN(iT) = r.s.H3PO4(iT) .* lim.C(iT) .* lim.P(iT) .* lim.NS(iT) ...
            .* (1 - lim.RS(iT)) .* (1/tau.fix) .* v.const.NPratio .* tdep.photo(iT); % mol N/yr
    end
% Prior to the advent of oxygenic photosynthesis, excess iron reduces N2 to
% NH3 during photo-ferrotrophy (anoxygenic fixation):
%                 472FeO + 106CO2 + 8N2 + H3PO4 + 1322H2O
%                                  -->
%             472Fe(OH)3 + (CH2O)_106(NH3)_16(H3PO4) + 484H2O
%    
% anoxygenic nitrogen fixation == 'fixation'
    if r.s.FeO(iT) > 0
        flux.fixation.FeO(iT) = r.s.FeO(iT) .* lim.C(iT) .* lim.P(iT) .* lim.NS(iT) ...
            .* (1 - lim.RS(iT)) .* (1/tau.fix) .* (1 - tdep.photo(iT)); % mol Fe/yr
    end
% both reactions assume immediate uptake of NH3 produced by fixation into
% organic matter (ie. not released into RN reservoir)
end

% Species produced and consumed during photo-ferrotrophic assimilation
sps  = {'CO2','H2O','RN','H3PO4'};                      
rats = [1,7,1/v.const.CNratio,1/v.const.CPratio]./4;        % ratios for each species per mole Fe used
for ips = 1:length(sps)
   flux.ferrotrophy.(sps{ips}) = flux.ferrotrophy.FeO .* rats(ips); 
end
% Species produced and consumed during oxygenic production of new fixed N 
fxlst = {'CO2','H2O','H3PO4','O2'}; 
raf    = [106 ,130, 1, 118]./16;                            % ratios for each species per mol newN used
for iaf = 1:length(fxlst)
   flux.fix.(fxlst{iaf}) = flux.fix.newN .* raf(iaf); 
end 
% Species produced and consumed during anoxygenic production of new fixed N 
afxlst = {'N2','H2O','CO2','H3PO4'}; 
rfx   = [8, 838, 106, 1]./472;                              % ratios for each species per mol Fe used
for ixf = 1:length(afxlst)
   flux.fixation.(afxlst{ixf}) = flux.fixation.FeO .* rfx(ixf); 
end

%--------------------------- Total Productivity ---------------------------
% The net reaction between direct oxygenic photosynthesis (assim), oxygenic
% fixation (fix), anoxygenic fixation (fixation), and photo-ferrotrophy.
%                                                       *** footnote 2 ***
flux.prod.CO2 = flux.assim + flux.fix.CO2...
    + flux.fixation.CO2 + flux.ferrotrophy.CO2;             % mol C/yr
flux.prod.H3PO4 = (flux.assim ./ v.const.CPratio) + flux.fix.H3PO4...
    + flux.fixation.H3PO4 + flux.ferrotrophy.H3PO4;         % mol P/yr
flux.prod.NH4 = vr.rassim.NH4.*(flux.assim./v.const.CNratio).*fRN;% mol N/yr; NH4 assimilation (for Gruber formalism)
flux.prod.NH3 = vr.rassim.NH3.*(flux.assim./v.const.CNratio).*fRN;% mol N/yr; NH3 assimilation (for canonical formalism)
flux.prod.RN  = flux.prod.NH4 + flux.prod.NH3;              % mol N/yr
flux.prod.N2  = flux.fixation.N2 + 0.5.*flux.fix.newN;      % mol N2/yr
flux.prod.HNO3 = (flux.assim ./ v.const.CNratio) .* (1-fRN);% mol N/yr
flux.prod.H2O = (flux.prod.HNO3.*vr.rassim.H2ONO)...
    + (flux.prod.RN.*vr.rassim.H2ONH) + flux.fix.H2O...
    + flux.fixation.H2O + flux.ferrotrophy.H2O;             % mol H2O/yr
flux.prod.FeO = flux.ferrotrophy.FeO + flux.fixation.FeO;   % mol Fe/yr
flux.prod.O2  = (flux.prod.HNO3.*vr.rassim.O2NO)...
    + (flux.prod.RN.*vr.rassim.O2NH) + flux.fix.O2;         % mol O2/yr; depending on which OM formalism is used, NH4 or NH3 flux will = 0

%------------------------ Terrestrial Productivity ------------------------
% A terrestrial biosphere evolves during the Phanerozoic, activating in 2 
% steps, first with fungal evolution (600-700 Ma), then with the evolution 
% of vascular land plants (400 Ma). This net flux is added/subtracted from
% atmospheric O2/CO2 reservoirs in the ODEs and org C is added to the 
% continent as CH2O (we assume the terrestrial biosphere is not strongly 
% nutrient limited).                                    *** footnote 3 ***
%                  CO2(a) + H2O(a) -hv-> CH2O(c) + O2(a)
% A small portion of organic carbon (10%) is not buried but instead
% remineralized by methanogens:
%                    10/100( CH2O --> 1/2 CH4 + 1/2 CO2 )
% We assume this flux has the same dependency on pCO2 as continental 
% weathering by vascular plants (1/2 exponent) since they are the source of
% this carbon, which also makes this somewhat responsive to temperature. 

flux.terrprod = 1.5e13 .* tdep.terr .* sqrt(r.a.CO2./v.atm.CO2pal);% mol O2/yr

%%                        Biogeochemical Reactions
%--------------------------------- Death ---------------------------------- 
% Your time is up, beasties!
flux.death.OC = zeros(refsize); 
for it = 1:length(t) 
    if r.s.LB(it) > 0
        flux.death.OC(it) = r.s.LB(it) ./ tau.death;        % mol C/yr
    end
end
np = {'ON','OP'}; 
for in = 1:length(np)
   flux.death.(np{in}) = flux.death.OC ./ v.const.(['C',np{in}(end),'ratio']); % mol/yr; org N and P generation by death  
end

% -------------- Nitrification and Biomass Remineralization ---------------
% These reactions all depend on dissolved oxygen concentrations and other
% environmental conditions, so we caculate them all in one convenient loop.
obx = {'s','d','n','z'}; 
for ii = 1:length(obx)
    rez = obx{ii}; 
    % initialize vectors
    flux.ammon.OC.(rez) = zeros(refsize); flux.denit.OC.(rez) = zeros(refsize); 
    flux.metha.OC.(rez) = zeros(refsize); flux.nitr.(rez) = zeros(refsize); 
    for iT = 1:length(t)
%---------------------------- Nitrification -------------------------------
% Requires oxygen and high pH (Gruber, 2008).
%              (NH4^+, NH3) + 2O2 --> HNO3 + (H^+) + H2O
        if r.(rez).O2(iT) > 0 && conc.(rez).NH4(iT) > 0
            flux.nitr.(rez)(iT) = r.(rez).RN(iT) .* lim.Onit.(rez)(iT) ... 
                .* (1-lim.Hnit.(rez)(iT)) .* (1./tau.nitr); % mol N/yr
        end
%----------------------- Biomass Remineralization -------------------------
% Pathway of choice depends on [O2] and [HNO3] conditions; flux magnitude
% depends on available OC. Anoxic pathways are less efficient and thus take
% longer (timescale multiplied by r_axrm, a constant). 
        % if [O2] > oxic limit, remineralize via ammonification
        if r.(rez).O2(iT) > 0 && r.(rez).OC(iT) > 0
            flux.ammon.OC.(rez)(iT) = (r.(rez).OC(iT) ./ tau.oxrm.(rez))...
                .* lim.Oaer.(rez)(iT);                      % mol C/yr
        end
        % given low [O2] but available HNO3, remin via heterophic denitrification 
        if r.(rez).OC(iT) > 0 && r.(rez).HNO3(iT) > 0
            flux.denit.OC.(rez)(iT) = (r.(rez).OC(iT) ./ (tau.oxrm.(rez) .* v.const.r_axrm))...
                .* (1-lim.Oaer.(rez)(iT)) .* lim.Nana.(rez)(iT);% mol C/yr
        end
        % given low [O2] and insufficient HNO3, remin via methanogenesis
        if r.(rez).OC(iT) > 0 
            flux.metha.OC.(rez)(iT) = (r.(rez).OC(iT) ./ (tau.oxrm.(rez) .* v.const.r_axrm))...
                .* (1-lim.Oaer.(rez)(iT)) .* (1-lim.Nana.(rez)(iT));% mol C/yr
        end
    end
end
spam = {'H3PO4','RN','H2O','O2'};                           % ammonification products
sphd = {'H3PO4','N2','RN','HNO3','H2O'};                    % denitrification products
spmt = {'H3PO4','CO2','CH4','RN','H2O'};                    % methanogenesis products

%---------------------- Aerobic Remineralization --------------------------
% Remineralization in open ocean is dependent on available oxygen (Quan + 
% Falkowski, 2009); we assume that organic nitrogen is
% transformed exclusively into its reduced form, ammonia/ammonium, and this 
% is then able to be nitrified if O2 is available.

%                           Ammonification
% Stoichiometry from Gruber (2008)
% C_106 H_175 O_42 N_16 P + 118O2 --> 106CO2 + 16NH4+ + HPO4 + 48H2O + 14OH-
% Canonical Stoichiometry (Fundamentals of Geobiology, Konhauser et al.): 
%   C(H2O)_106 (NH3)_16 H3PO4 + 106O2 --> 106CO2 + 16NH3 + H3PO4 + 106H2O
for id = 1:length(obx)
    rs = obx{id}; 
    for ip = 1:length(spam)
        flux.ammon.(spam{ip}).(rs) = flux.ammon.OC.(rs) .* vr.rammon.(spam{ip});% mol/yr
    end

%----------------------- Anaerobic Remineralization -----------------------
% In anoxic conditions all labile OM is remineralized without oxygen by two
% pathways (depending on HNO3 availability).

%                    Heterotrophic denitrification
% Stoichiometry from Gruber (2008):
% C_106 H_175 O_42 N_16 P + 104HNO3 --> 106CO2 + 60N2 + H3PO4 + 138H2O
% Canonical Stoichiometry (Fundamentals of Geobiology, Konhauser et al.):
% C(H2O)_106 (NH3)_16 H3PO4 + (424/5)HNO3 --> 106CO2 + (212/5)N2 + 16NH3 + H3PO4 + (742/5)H2O
    for is = 1:length(sphd)
        flux.denit.(sphd{is}).(rs) = flux.denit.OC.(rs) .* vr.rdenit.(sphd{is});% mol/yr
    end
    
%                           Methanogenesis  
% Stoichiometry from Gruber (2008):
% C_106 H_175 O_42 N_16 P + 59H2O --> 47CO2 + 59CH4 + 16NH3 + H3PO4 + 3H2O
% Canonical Stoichiometry (Fundamentals of Geobiology, Konhauser et al.):
%      C(H2O)_106 (NH3)_16 H3PO4 --> 53CO2 + 53CH4 + 16NH3 + H3PO4 

    for iv = 1:length(spmt)
        flux.metha.(spmt{iv}).(rs) = flux.metha.OC.(rs) .* vr.rmeth.(spmt{iv});% mol/yr
    end

% ---------------------------- Methanotrophy ------------------------------
% Using available oxygen, methanotrophs consume the methane produced in
% anoxic environments and convert it to CO2 (Goldblatt et al. 2006)
%                        CH4 + 2O2 --> CO2 + 2H2O
    flux.mtrophy.(rs) = zeros(refsize);
    for iT = 1:length(t)
        if r.(rs).CH4(iT) > 0 && r.(rs).O2(iT) > 0
            flux.mtrophy.(rs)(iT) = (r.(rs).CH4(iT) .* lim.Mtro.(rs)(iT)) ./ tau.mtrophy;% mol C/yr
        end
    end
end


%%                           Ocean Geochemistry
% ------------------------ Abiotic Iron oxidation -------------------------   
% Reduced iron is oxidized by available oxygen (disregarding pH here since 
% we assumed a neutral Fe(II) species). This flux follows the same 
% limitation as photoferrotrophy (lim.Phx) but in the opposite direction. 
% We use the following stoichiometry:
%              4FeO + O2 + 11H2O --hv--> 4Fe(OH)3 + 5H2O
% Rate constant Kphotox is from Braterman et al. (1986; 150-750 mg/cm2 yr)
flux.photox.FeO = zeros(refsize); 
for it = 1:length(t)
    if r.s.O2(it) > 0 && r.s.FeO(it) > 0
        flux.photox.FeO(it) = r.s.FeO(it) .* lim.Phx(it) .* (1/tau.photox);% mol Fe/yr
    end
end
sps = {'O2','H2O'}; ratxs = [1, 6]./4;                      % ratios of reactants/products
for ixs = 1:length(sps)
   flux.photox.(sps{ixs}) = flux.photox.FeO .* ratxs(ixs);  
end

%---------------------- Precipitation + Dissolution -----------------------
% Based on DIC/TA conditions in an ocean box, calcium carbonate is either
% precipitated (forward reaction) or dissolved (reverse reaction):
%                 Ca^2+ + 2HCO3^- <--> CaCO3 + CO2 + H2O
% Calcite/aragonite saturation (Omega - calculated in Flux_Spec) determines
% if precipitation or dissolution occurs. 
for iv = 1:length(obx)
    flux.precip.(obx{iv}) = zeros(refsize); flux.diss.(obx{iv}) = zeros(refsize); 
    for ip = 1:length(r.s.CaCO3)
        if conc.Oa.(obx{iv})(ip) >= 1  % saturation state
        % if saturation is reached, CaCO3 precipitates and doesn't dissolve
            flux.precip.(obx{iv})(ip) = ((conc.Oa.(obx{iv})(ip) - 1)^v.const.n_pre) *...
               v.oc.m.(obx{iv}) * conc.(obx{iv}).CO3(ip) ./ tau.precip;
            flux.diss.(obx{iv})(ip) = 0;
        elseif (conc.Oa.(obx{iv})(ip) < 1) && (conc.Oa.(obx{iv})(ip) > 0) && (r.(obx{iv}).CaCO3(ip) > 0)
        % if saturation is too low and CaCO3 is present, it dissolves
            flux.diss.(obx{iv})(ip) = (1 - conc.Oa.(obx{iv})(ip))^v.const.n_diss *...
               r.(obx{iv}).CaCO3(ip) ./ tau.diss.(obx{iv}); 
            flux.precip.(obx{iv})(ip) = 0;
        elseif (r.(obx{iv}).CaCO3(ip) <= 0) || (conc.Oa.(obx{iv})(ip) <= 0)
            flux.precip.(obx{iv})(ip) = 0; flux.diss.(obx{iv})(ip) = 0;
        end
    end
end

%------------------------- Reverse Weathering -----------------------------
% This "back reaction" of cations and alkalinity generates authigenic clays
% and acidity (Isson + Planavsky, 2018 eq. 3).          *** footnote 4 ***
% Stoichiometrically, this is treated:
%            Ca^2+ + SiO2 + 2HCO3^- --> CaSiO3 + 2CO2 + H2O
sedbox = {'n','z'};
for isb = 1:length(sedbox)
    bx = sedbox{isb}; 
    flux.revweather.(bx) = zeros(refsize); 
    krw.(bx) = 1.01e-19 .* conc.(bx).pH.^22.4;              % /yr; rate constant based on pH
    % Dissolved [SiO2] is kept near the upper estimate for modern (<0.1 mM)
    % Converted to units of mol eq/yr aka mol C/yr because we don't track Si
    for ita = 1:length(t)
        if r.(bx).TA(ita) > 0 
            flux.revweather.(bx)(ita) = 2.*krw.(bx)(ita).*inp.dSi.*v.oc.vol.(bx);% mol TA eq/yr
        end
    end
end

%------------------------ Phosphate Adsorption ----------------------------
% Iron oxide minerals absorb dissolved phosphate in certain conditions,
% forming iron oxide-phosphates: 
%                   Fe(OH)3 + H3PO4 --> FePO4 + 3H2O
% This requires a pH < 9 (Bjerrum + Canfield, 2022; Ajmal et al. 2018) so
% that phosphate speciates from H3PO4 -> PO4^3-, and is rate controlled by 
% the concentrations of P, Fe(OH)3 and SiO2, and a distribution coefficient
% K_D; Planavsky et al. 2017). For simplicity, we ignore the silica 
% constraints here. From Bjerrum + Canfield:
%                      F = K_ads .* [PO4] .* xF_iron
% xF_iron is the fraction of iron influx as oxides; we instead use Fe(OH)3
% reservoir over a sorption reaction timescale. 
for iss = 1:length(obx)
    flux.sorb.(obx{iss}) = (r.(obx{iss}).FeOH3./tau.sorb)...
        .* (1-lim.Hsb.(obx{iss})) .* (c.(obx{iss}).H3PO4 .* v.KdP);% mol P/yr
end

%------------------------- Seafloor Weathering ----------------------------
% Similar to subaerial silicate weathering, this is a direct control on
% atmospheric CO2 (Sleep + Zahnle, 2001; Mills et al. 2014). As pCO2
% increases, carbonate minerals may be preciptated in the porespaces of sea
% seafloor basalts:
%                     CaSiO3 + CO2 --> CaCO3 + SiO2
% CO2 control is generally related by power law to the relative size of the 
% atmospheric CO2 reservoir (with a weak control, a = 0.23; Mills et al 
% 2014). We assume that this takes on the same form as our other 
% CO2-dependant modifiers, but with a normalization of to spreading rate
% (fluid flow) in the modern context, and DIC dependence. Estimated modern
% rate from Mills et al. (2014) and Gillis + Coogan (2011). 
flux.sfw = 5e12 .* ws.sfw;                                  % mol C/yr

%--------------------- Hydrothermal Sequestration -------------------------
% A small portion of dissolved NH4 replaces K in alkali feldspars during
% hydrothermal alteration, in a reverse of silicate-N weathering:
% 4NH4+(aq) + 2Al2Si2O5(OH)4 + 8SiO2 + 4HCO3^-(aq) --> 4NH4AlSi3O8(s) + 4CO2(aq) + 6H2O
% simplified to reduced untracked species:
%        2NH4^+(aq) + SiO2 + 2HCO3^- -> (NH4)2SiO3 + 2CO2(aq) + H2O 
% Similar to seafloor weathering, we assume that this reaction is
% sensitive to CO2 and T conditions (which control SiO2 release) and 
% removes NH4 from the reactive sediments, incorporating it into oceanic 
% crust. The timescale for this reaction is a function of hydrothermal flow
% estimates for the modern day.

hyd = 1.6e6.*v.conv.spyr;                                   % m3/yr; modern hydrothermal volumetric flow (from 1.6 sverdrups; Johnson + Goldblatt, 2017)
flux.hyd = c.z.NH4 .* hyd .* ws.sfw;                        % mol N/yr

%%                            Mass Movement
% ----------------------- Export + Sedimentation --------------------------
% Export is from surface ocean to deep ocean ONLY. While shelf area 
% fraction (f_Ashelf) is small, the amount of productivity in shelf vs open
% ocean (f_Pshelf) is substantial; therefore most OM and CaCO3 is
% sedimented into the neritic sediments rather than exported to the deep 
% ocean. All deep ocean solids are sedimented in deep sediments. OM, CaCO3, 
% and iron oxides (BIFs) sink out along different timescales. 
sps = {'OC','ON','OP','CaCO3','FeOH3','FePO4'}; 
for is = 1:length(sps)
    fshelf = f_Pshelf; 
    switch sps{is}
        case {'FeOH3','FePO4'} % use shelf area and bifsink timescale
            tau_S = tau.bifsink.s;  tau_D = tau.bifsink.d; fshelf = v.f.Ashelf; 
        case 'CaCO3' % use carbsink timescale
            tau_S = tau.carbsink.s; tau_D = tau.carbsink.d; 
        otherwise
            tau_S = tau.sink.s;     tau_D = tau.sink.d; 
    end
    flux.export.(sps{is}) = ((1-fshelf) .* r.s.(sps{is})) ./ tau_S; % mol/yr; export to deep ocean
    flux.sed.(sps{is}).n = (fshelf .* r.s.(sps{is})) ./ tau_S;% mol/yr; sedimentation into shallow seds
    flux.sed.(sps{is}).z = r.d.(sps{is}) ./ tau_D;          % mol/yr; sedimentation into deep seds
end

%------------------------- Sedimentary Burial -----------------------------
% Reactive sediments (n,z) are removed from interactions with the ocean, 
% thus becoming unreactive sediments (u). Sediment residence time (burial 
% timescale) is proportional to erosional flux from continents (xeros), 
% which we assume here is proportional to continental size (ie. larger 
% continents = more erosional surface = higher sediment transport to ocean,
% faster burial). NOTE: burial of FePO4 takes the form: 
%              FePO4 + 3HCO3^- --> Fe(OH)3 + PO4^3- + 3CO2
% requiring -3TA per mol FePO4

for izn = 1:length(sedbox)
    sbx = sedbox{izn};  
    dsp = {'OC','CaCO3','OP','ON','FeOH3','FePO4'};
    for is = 1:length(dsp)
        flux.burial.(dsp{is}).(sbx) = r.(sbx).(dsp{is}) ./ (tau.sink.(sbx) .* xeros);% mol/yr
    end
    % ------------------ P Sequestration in Carbonate ---------------------
    % Some dissolved phosphate gets bound to calcium carbonate as it 
    % precipitates. Stoichiometrically this is treated as switching with 
    % the carbonate ion rather than making Ca5(PO4)3F:
    %             2H3PO4 + 3CaCO3 --> Ca3(PO4)2 + 3CO2 + 3H2O
    % This flux depends on the burial flux of CaCO3 and relative [H3PO4]
    % with a constant fraction for how much apatite is found in general
    % carbonate rock                                    *** footnote 5 *** 
    xP.(sbx) = c.(sbx).H3PO4 ./ (v.oc.dP.(sbx) .* v.oc.rho);% [H3PO4] relative to modern estimates (concentrations in mol/m3)
    flux.burial.PO4.(sbx) = (2/3) .* flux.burial.CaCO3.(sbx)...
        .* xP.(sbx) .* v.f.apatite;                         % mol P/yr
end

% ---------------------------- P Scavenging -------------------------------
% Laminated (anoxic) shales tend to have higher organic C:P ratios,
% suggestive of P scavenging in phosphorus-limited conditions by
% remineralizing bacteria (Ingall et al. 1993, VanCapellen + Ingall, 1994).
%     (CH2O)_106(NH3)_16(H3PO4)(s) -> 106CH2O(s) + 16NH3(s) + H3PO4(aq)
% We assume that this happens at a rate proportional to the degree of 
% anoxicity and to bio P scarcity. This flux can produce buried OM C:P 
% nearly 1 mag greater than RR (Alcott et al. find ~ 1100:1, Van Cappellen
% + Ingall 1994 find 3900:1 max) - we max total scavenging to a fixed
% fraction (f_scavP) of organic P. This reaction occurs at a rate
% controlled by reactive sediment residence. 
for ig = 1:length(sedbox)
    fscavP.(sedbox{ig}) = zeros(refsize); 
    for io = 1:length(c.(sedbox{ig}).O2)
        if (r.(sedbox{ig}).OP(io) > 0) 
            fscavP.(sedbox{ig})(io) = v.f.scavP .* (1 - lim.P(io)) .* (1 - lim.Anox.(sedbox{ig})(io));% fraction of P removed
        else
            fscavP.(sedbox{ig})(io) = 0;  
        end
    end
    flux.scav.(sedbox{ig}) = fscavP.(sedbox{ig}) .* (r.(sedbox{ig}).OP ./ (tau.sink.(sedbox{ig}) .* xeros));
end

%%                          Downgoing Material
% ------------ Sediment/Slab Accretion, Volcanism, Subduction -------------
% Sediments and oceanic crust are exported to the continents, atmosphere,
% or mantle via accretion, volcanism, or subduction fluxes respectively.
% All of these operate along the same timescale; the canonical lifetime of 
% the oceanic crust (tau_subd).                         *** footnote 6 ***
% NOTE: P species are not melted, but are instead transferred to the 
% silicate reservoirs on the continent.  

% We use fixed fractions which sum to 1 (ie. f_acc + f_melt + f_sub == 1)
% for all pathways, with the exception of some tests for variable accretion
% fraction.                                             *** footnote 7 ***
f_acc.u = 2.*f_acc.o; % potentially variable; f_acc for ocean crust defined above 
% Variable accretion fraction only effects suduction fraction, melt is kept
% constant.
uobx = {'u','o'}; 
for ib = 1:length(uobx)
    rs = uobx{ib}; 
    sps = fieldnames(r.(rs)); 
    f_sub.(rs) = 1 - (f_acc.(rs) + v.f.melt);               % subducted fraction based on accreted fraction
    for is = 1:length(sps)
        spe = sps{is}; nm = spe; 
        f_volc = v.f.volc.(nm);                             % volcanic fraction specific to this species (portion of melt)
        f_cryst = v.f.melt - f_volc;                        % recrystallized fraction (other portion of melt)
        switch spe
            case 'CaCO3' % change in name denotes source reservoir (carb == sediments, CaCO3 == slab)
                if strcmp(rs,'u') 
                    nm = 'carb';
                end  
        end
        flux.acc.(nm)       = (r.(rs).(spe) ./ tau.subd) .* f_acc.(rs);% mol/yr; accretion flux
        flux.subduct.(nm)   = (r.(rs).(spe) ./ tau.subd) .* f_sub.(rs);% mol/yr; subduction flux
        if f_volc > 0 % if these fractions are 0, don't even make a flux 
            flux.volc.(nm)  = (r.(rs).(spe) ./ tau.subd) .* f_volc;% mol/yr; volcanism flux
        end
        if f_cryst > 0
            flux.cryst.(nm) = (r.(rs).(spe) ./ tau.subd) .* f_cryst;% mol/yr; recrystallization flux
        end
    end
end
% organic carbon is volatilized as CO2 and CH4 in equal proportions
%                    CH2O -heat-> 1/2 CO2 + 1/2 CH4
flux.volc.CO2 = 0.5.*flux.volc.OC ;                         % CHONP -> CO2
flux.volc.CH4 = flux.volc.CO2;                              % CHONP -> CH4

%---------------------------- Metamorphism --------------------------------
% Continental materials exposed to heat undergo contact metamorphism and
% release volatiles to the atmosphere. We assume this material has a longer
% residence time in the continents than in the arc wedge (ie. tau_meta > 
% tau_sub), and is also related to spreading (faster spreading being 
% related to continental uplift, more accretion events, etc). Same
% stoichiometries as volcanism for the following species. 
% NOTE: metamorphosed P species are transferred to carbonates 
% (Ruttenberg, 2003).

spc = {'OC','ON','OP','CaCO3','NH4'};
for il = 1:length(spc)
    flux.meta.(spc{il}) = (r.c.(spc{il}) ./ tau.meta);      % mol/yr
end   
flux.meta.CO2 = 0.5.*flux.meta.OC ;                         % CHONP -> CO2
flux.meta.CH4 = flux.meta.CO2;                              % CHONP -> CH4

end

%% Subfunction : Determine assumed organic matter stoichiometry

function [vr] = SetStoichiometricRatios(STOINM)

% make the stoichiometric balances for different production/remin pathways
% adjustable according to the organic matter stoichiometry setting and
% calculated P sensititivty
switch lower(STOINM)
    case 'gruber'

    % ASSUMING GRUBER (2008) OM STOICHIOMETRY : C_106 H_175 O_42 N_16 P
    rx.assim.O2NO    = 150/16;           % O2 produced/NO3; assimilation O2 ratio if NO3 utilized 
    rx.assim.O2NH    = 118/16;           % O2 produced/NH4; assimilation O2 ratio if NH4 utilized
    rx.assim.H2ONO   = 78/16;            % H2O/NO3; assimilation H2O uptake ratio if NO3 utilized
    rx.assim.H2ONH   = 48/16;            % H2O/NH4; assimilation H2O uptake ratio if NH4 utilized
    rx.assim.HNO     = 18/16;            % H+/NO3; assimilation H uptake when NO3 used
    rx.assim.OHNH    = 14/16;            % OH-/NH4; assimilation OH uptake when NH4 used
    rx.assim.NH4     = 1;                % this assimilation formulation uses NH4+, not NH3
    rx.assim.NH3     = 0;                % ^^
    rx.ammon.H3PO4   = 1/106;            % H3PO4/total C; ammonification PO4 release ratio
    rx.ammon.RN      = 16/106;           % NH4 produced/total C; ammonification NH4 prod ratio
    rx.ammon.H2O     = 48/106;           % H2O produced/total C; ammonification H2O prod ratio
    rx.ammon.OH      = 14/106;           % OH produced/total C; ammonification OH prod ratio
    rx.ammon.O2      = 118/106;          % O2/total C; ammonification O2 uptake ratio 
    rx.denit.H3PO4   = 1/106;            % HPO4/total C; hetero-denit PO4 release ratio
    rx.denit.N2      = 60/106;           % N2/total C; hetero-denit N2 prod ratio
    rx.denit.RN      = 0/106;            % NH3/total C; hetero-denit NH3 prod ratio (ie. not produced!)
    rx.denit.HNO3    = 104/106;          % NO3/total C; hetero-denit NO3 uptake ratio
    rx.denit.H2O     = 138/106;          % H2O produced/total C; hetero-denit H2O prod ratio
    rx.denit.H       = 104/106;          % H uptake/total C; hetero-denit H uptake ratio
    rx.meth.H3PO4    = 1/106;            % HPO4/total C; methanogenesis PO4 release ratio
    rx.meth.CO2      = 47/106;           % CO2/total C; methanogenesis CO2 prod ratio
    rx.meth.CH4      = 59/106;           % CH4/total C; methanogenesis CH4 prod ratio
    rx.meth.RN       = 16/106;           % NH3/total C; methanogenesis NH3 production
    rx.meth.H2O      = 56/106;           % H2O/total C; methanogenesis H2O uptake ratio (59 mol in -> 3 mol out)

    case 'canon'
    % ASSUMING "CANONICAL" OM STOICHIOMETRY : C_106 H_263 O_110 N_16 P
    % from Fundamentals of Geobiology Konhauser et al. + Berman-Frank et al. (2008)
    rx.assim.O2NO    = 138/16;           % O2 produced/NO3; assimilation O2 ratio if NO3 utilized 
    rx.assim.O2NH    = 106/16;           % O2 produced/NH3; assimilation O2 ratio if NH3 utilized
    rx.assim.H2ONO   = 122/16;           % H2O/NO3; assimilation H2O uptake ratio if NO3 utilized
    rx.assim.H2ONH   = 106/16;           % H2O/NH3; assimilation H2O uptake ratio if NH3 utilized
    rx.assim.HNO     = 16/16;            % H+/NO3; assimilation H uptake when NO3 used
    rx.assim.OHNH    = 0/16;             % OH-/NH3; assimilation OH uptake when NH3 used
    rx.assim.NH4     = 0;                % this assimilation formulation uses NH3, not NH4+
    rx.assim.NH3     = 1;                % ^^
    rx.ammon.H3PO4   = 1/106;            % H3PO4/total C; ammonification PO4 release ratio
    rx.ammon.RN      = 16/106;           % NH3 produced/total C; ammonification NH3 prod ratio
    rx.ammon.H2O     = 106/106;          % H2O produced/total C; ammonification H2O prod ratio
    rx.ammon.OH      = 0/106;            % OH produced/total C; ammonification OH prod ratio
    rx.ammon.O2      = 106/106;          % O2/total C; ammonification O2 uptake ratio 
    rx.denit.H3PO4   = 1/106;            % H3PO4/total C; hetero-denit PO4 release ratio
    rx.denit.N2      = 42.4/106;         % N2/total C; hetero-denit N2 prod ratio
    rx.denit.RN      = 16/106;           % NH3/total C; hetero-denit NH3 prod ratio
    rx.denit.HNO3    = 84.8/106;         % NO3/total C; hetero-denit NO3 uptake ratio
    rx.denit.H2O     = 148.4/106;        % H2O produced/total C; hetero-denit H2O prod ratio
    rx.denit.H       = 84.8/106;         % H uptake/total C; hetero-denit H uptake ratio
    rx.meth.H3PO4    = 1/106;            % H3PO4/total C; methanogenesis PO4 release ratio
    rx.meth.CO2      = 53/106;           % CO2/total C; methanogenesis CO2 prod ratio
    rx.meth.CH4      = 53/106;           % CH4/total C; methanogenesis CH4 prod ratio
    rx.meth.RN       = 16/106;           % NH3/total C; methanogenesis NH3 prod ratio
    rx.meth.H2O      = 0/106;            % H2O/total C; methanogenesis H2O uptake ratio


end

vr.rassim = rx.assim; 
vr.rmeth  = rx.meth;
vr.rdenit = rx.denit;
vr.rammon = rx.ammon; 
end

%% Subfunction: check real time in model run

% inputs are model time vector (t) and the input structure (inp). Outputs a
% new time vector with "true" model time, if this is a reset run.

function nt = CheckTheTime(t,inp)
    nt = t;                                                 % baseline assumption that time is valid
    if isfield(inp,'treset')
        nt = t + inp.treset;                                % add the reset time to the time vector we're evaluating
    end

end


%% Subfunction: evolve parameter of choice
% Technically, in the nominal run this is only used to linearly-decrease
% mantle reductant influx, but in theory can be used to linearly change any
% tunable parameter of choice!

function [evp] = EvolveParam(time,endt,pointA,pointB,refsz)
% linearly evolve a parameter between two endmembers (inp.A = starting point,
% inp.B = end point) across a certain portion of the model run (endt)
pert = zeros(refsz);                  % initialize perturbation timing
for it = 1:length(time)
     % the amount of parameter change is proportional to time between point A and B
    if time(it) <= endt 
        pert(it) = time(it)./endt; 
    else % after end time, the parameter is kept at the point B value
        pert(it) = 1; 
    end
end
evp = pointA + ((pointB - pointA).*pert); 

end

%% ---------------------------- FOOTNOTES ---------------------------------
%                               *** 1 *** 
% The nitrogenase enzyme likely evolved early (Stuken et al. 2015; Zerkle 
% et al. 2006) and primarily uses molybdenum as its most efficient metal
% cofactor (can use vanadium, iron, mixes of all three). Prior to oxidative
% weathering, Mo was less abundant in the ocean and likely limited the
% spread of this particular pathway. After oxygen rose, fixing may have
% been Mo limited because of higher prevalence of oxic sinks (Thoby et al.
% 2019). Conversely Zerkle et al. (2006) argue that rise in O2 led to 
% higher [Mo] and therefore higher productivity, regardless of nitrogenase
% enzyme deactivation in presence of O2 - temporal, spatial, or reductive
% techniques to separate nitrogenase from environmental O2 seem to have
% evolved early in cyanobacteria. Nitrogenase appears to have first evolved
% in methanotrophic bacteria/archea, and preursor (and much less efficient) 
% nitrogen fixation methods appear to have evolved early in many different
% types of anaerobic organisms (Boyd + Peters, 2013).
% Given the multitude of fixation methods throughout the biosphere and the
% seeming ubiquity of this function in primary producers, we assume that 
% fixation has always been a component of the biosphere. In addition, we
% make no assumption regarding sensitivity to oxygen at low or high
% concentrations, again since there are so many different nitrogenase
% alternatives with different oxygen sensitivies, many methods by which
% organisms overcome nitrogenase deactivation (Stuken et al 2015), and so 
% many overlapping pathways for fixation that don't use nitrogenase or 
% associate metal cofactors Boyd + Peters, 2013). In sum, nitrogen fixation
% seems to be such an important biological function that if one pathway in
% the ecosystem becomes inhibited by oxygen, it's more than likely replaced 
% by an althernative fixation mechanism somewhere else in the biosphere.

%                               *** 2 ***
% The total production uptake fluxes for RN species do NOT include 
% ferrotrophy because of knock-on effects for H2O uptake calculation w/r/to 
% photosynthetic assimilation of HNO3 vs. NH3/NH4. Photo-ferrotrophic
% uptake of RN is done in the ODEs instead.

%                               *** 3 ***
% Organic carbon burial on land is somewhere between 2.5 - 4.5e12 mol C/yr 
% (Lenton et al, 2018) and 350e12 g C/yr (~3e13 mol C/yr; Berner, 2009) -
% we assume an intermediate value between these two as our fixed flux. 
% Based on IPCC estimates, this flux should broadly equal modern day OC 
% burial in the ocean. 

%                               *** 4 ***
% Reverse weathering is further controlled by cation source and the 
% dissolved concentration of silica, which in particular has changee over
% geologic time as siliceous organism (diatoms) evolve; this development 
% resulted in a decrease in the efficiency of this buffer (Isson + 
% Planavsky, 2018). For the sake of stoichiometric closure, this is treated 
% as a direct reversal of the subaerial silicate weathering flux (ie. in
% I+P 2018 equation 3, cation X^2+ = Ca^2+) and we don't assume cation
% source (ie. hydrothermal or from detrital clays) because we don't 
% explicitly track Ca^2+. For the sake of simplicity, this flux is not 
% dependent on calculated porewater [SiO2]; we have simplified the pH 
% dependencies from I+P 2018 because the Si0 dependency is much weaker 
% than the rate constant (k_rw) pH dependency anyway. An additional forcing
% in future versions of this model could be a change in [SiO2] after the
% evolution of diatoms.

%                               *** 5 *** 
% Kraal et al. (2017) estimate 10-15% Ca-bound P in detrital sediments;  
% organic and carbonate bound-P together make up ~ 30% of the total deep 
% sediment P reservoir. Knudsen + Gunter (2002) measured the prevalence of 
% apatite in sedimentary carbonates (phosphorites) in Ordovician limestone
% (5-6 wt%). We use a conservative estimate of 1% as our fraction of P 
% associated with CaCO3 precipitation nand assume this burial flux grows if 
% the phosphorus reservoir increases over modern concentrations 
% (which is taken from WOCE atlas dataset from UC Davis). This fraction
% results in a modern inorganic P burial flux (1-5e11 mol/yr) under modern
% [H3PO4] and sedimentary carbonate burial conditions. 

%                               *** 6 ***
% The stoichiometries for downgoing material are detailed, with untracked
% species denoted by (nt) comments- these are important to add to the mass
% conservation reservoirs for O and H in the ODEs. 
% 
% Stoichiometries for volcanism/metamorphism:
% NH4:            1/2(NH4)2SiO3 -heat-> NH3 + 1/2H2O(a) + 1/2SiO2 (nt)                         
% OC,ON,OP:
%  CH2O_106 NH3_16 H3PO4 -heat-> 53CO2(a) + 53CH4(a) + 16NH3(a) + PO4(c) + 3H (nt)
% CP:             1/2Ca3(PO4)2 -->  PO4 + 3/2Ca (nt)
%  OP + CP extra step: PO4 + 1/2Fe2SiO4 --> FePO4 +  1/2SiO4 (nt)
%  ON recrystallization: 
%                 2NH3 + SiO2 + H2O(u) --> (NH4)2SiO3(c) 
% carb,CaCO3:     CaCO3 + SiO2 -heat-> CaSiO3 + CO2(a)
% 
% Stoichiometries for subduction:
% CP,FePO4:       1/2Ca3(PO4)2, PO4 --> P (Ca ignored)
% OC,ON,OP: 
% CH2O_106 NH3_16 H3PO4 --> 106C + 16N + 1P + 263H (nt) + 212O (nt)
% NH4:            NH4 --> N
% FeOH3, Fe2SiO4: (FeOH3, Fe2SiO4) --> (1,2)Fe + (3H + 3O, 4O) (nt)
% 
% No change in stoichiometries for accretion, transfers directly to a 1:1 
% stoichiometric continental reservoir. 

%                               *** 7 ***
% We assume sediments are 2x as likely to be accreted (potentially up to
% 70% of sediments forming accretionary prism in modern subduction zones;
% von Heune + Scholl, 1991; Wong et al. 2019). This fraction is not well 
% constrained for the distant past, therefore we take a very conservative 
% estimate of up to 10% sediments accreted across all time. 
% In some tests, we may evolve f_acc (and thus the other fractions to == 1 
% in sum) relative to continental area. 
% 
% Melt fraction is always fixed at same value, but for different species
% the resulting pathways are split differently between recrystallization
% onto the continent and volcanic outgassing (f_crys + f_volc = f_melt). 
% 
% Modern volatile subduction fractions are 19% for N (from Mallik et al. 
% 2018), and 0.5% of total slab C (based on the nearly complete 
% volatilization found by Li et al. 2020). A cooler slab and mantle both 
% contribute to enhanced retention of volatiles, therefore one could relate
% the subducted fraction (f_sub) to mantle temperature (Johnson + Goldblatt,
% 2017). Here, we keep subduction fractions fixed, and opposite to
% accretion (ie. sediments are half as likely to be subducted as the
% oceanic crust). 

