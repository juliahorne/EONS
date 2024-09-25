%% ============================= EONS Model =============================== 
% Julia Horne, 2018

% Includes constant values for calculations, conversions, fractions, 
% earth/atm/ocean measurable properties, and color assignments for 
% plotting. Also reads and compiles the output from Byrne 2014's 
% line-by-line transfer model and sets up a constant structure
% (ie. bb) to reference this dataset for interpolation in RFInterp
% function. Don't worry! It all comes out as a single structure (v). 

set(groot,'defaultLineLineWidth',2); % for plots

% -------------------------- model constants ------------------------------

const.mm.CaSiO3= 0.116;            % kg/mol; molar mass of wollastonite
const.mm.SiO3  = 0.076;            % kg/mol; molar mass of silicate
const.mm.sil   = const.mm.CaSiO3;  % kg/mol; average molar mass of relative Ca ~4.25 wt% + Mg 2.5 wt% silicate crust
const.mm.N     = 0.014;            % kg/mol; molar mass of Nitrogen  
const.mm.P     = 0.031;            % kg/mol; molar mass of phosphorus
const.mm.Fe    = 0.05585;          % kg/mol; molar mass of iron 
const.mm.H2O   = 0.018;            % kg/mol; molar mass of water
const.mm.NH4   = 0.018;            % kg/mol; molar mass of ammonium
const.mm.NH3   = 0.017;            % kg/mol; molar mass of ammonia
const.mm.NO3   = 0.062;            % kg/mol; molar mass of nitrate
const.mm.PO4   = 0.095;            % kg/mol; molar mass of phosphate
const.mm.CO2   = 0.044;            % kg/mol; molar mass of carbon dioxide
const.mm.CO3   = 0.060;            % kg/mol; molar mass of carbonate
const.mm.CaCO3 = 0.101;            % kg/mol; molar mass of Calcite
const.mm.C     = 0.012;            % kg/mol; molar mass of carbon
const.mm.CH4   = 0.016;            % kg/mol; molar mass of methane
const.mm.O2    = 0.032;            % kg/mol; molar mass of O2 gas
const.mm.O     = 0.016;            % kg/mol; molar mass of oxygen
const.mm.FeOH3 = 0.107;            % kg/mol; molar mass of iron oxide
const.g        = 9.81;             % m/s2; gravity
const.R        = 8.3145;           % m^3Pa/Kmol; ideal gas law constant (aka J/Kmol)
const.Rhat     = 83.14510;         % cm^3bar/Kmol; ideal gas law constant for pressurized equilibrium constants (A.11, Zeebe + Wolf-Gladrow)
const.bol      = 5.67e-8;          % W/m^2K^4; Stefan-Boltzmann constant
const.NPratio  = 16/1;             % ratio of N to P according to Redfield Ratio OM = C_106 H_263 O_110 N_16 P
const.CNratio  = 106/16;           % ratio of C to N according to Redfield Ratio OM = C_106 H_263 O_110 N_16 P 
const.CPratio  = 106/1;            % ratio of C to P according to Redfield Ratio OM = C_106 H_263 O_110 N_16 P 
const.COratio  = 106/110;          % ratio of C to O according to Redfield Ratio OM = C_106 H_263 O_110 N_16 P
const.CHratio  = 106/263;          % ratio of C to H according to Redfield Ratio OM = C_106 H_263 O_110 N_16 P
const.OFeratio = 18/12;            % ratio of O to Fe in iron oxidation reaction (ie. O:Fe = 18:12)
const.CS       = 0.8;              % W/m^2/K; climate sensitivity parameter (Byrne + Goldblatt, 2014)
const.tauT     = 50;               % yr; ocean-atmosphere equilibration timescale
const.rho.gran = 2700;             % kg/m3; density of granite
const.rho.basa = 3000;             % kg/m3; density of basalt
const.r_axrm   = 10;               % ratio of anoxic:oxic remineralization rates (est 5-50:1) 
const.n_diss   = 1.05;%2.4;        % dissolution order (tuned from LOSCAR quoted value- on right)
const.n_pre    = 1.25;%4.5;        % precipitation order (see ^^)
    
% ------------------------- conversion values -----------------------------

conv.hryr      = 365*24;           % hours/yr; duration of 1 year in hours
conv.spyr      = 60*60*24*365;     % s/yr; years to seconds conversion factor
conv.L2m3      = 1000;             % L/m^3; liters per cubic meter volume conversion factor
conv.Kelv      = 273.15;           % degrees K @ 0C
conv.L2kg      = 1027/1000;        % kg/L; L to kg conversion factor
conv.diff      = 60*60*24*365/100^2;% m2/yr; cm^2 s^-1 conversion for diffusion
conv.Jcal      = 4.182;            % J/cal; convert Joules to calories
conv.cm3m3     = 100^3;            % cm^3/m^3; convert cm^3 to m^3
conv.cm2m2     = 100^2;            % cm^2/m^2; convert cm^2 to m^2
conv.barPa     = 1/100000;         % bar/pascal; convert pressure pascal to bar (1 pascal = 1 N/m2 = 1 kg/ms2)
conv.depthbar  = conv.barPa*const.g*1027;% m/bar; convert depth to bar
conv.depthatm  = conv.depthbar.*1.01325; % m/atm; convert depth to atm 
conv.avo       = 6.02e23;          % molecules/mol; avogadro's number 
conv.ppm       = 1e6./1.773e20;    % ppm; parts per million in atm
conv.ppb       = 1e9./1.773e20;    % ppb; parts per billion in atm
conv.atmunits  = (conv.spyr*conv.cm2m2*5.101e14)/conv.avo;% convert stupid molecules/cm2/yr to mol/yr

% ----------------------------- atmosphere --------------------------------
ea.sa          = 5.101e14;         % m^2; earth surface area
ea.rad         = 6378e3;           % m; equatorial radius of earth
atm.mm         = 0.03;             % kg/mol; approx. molar mass of atmosphere based on composition
atm.mol        = 1.773e20;         % number of moles in atmosphere
atm.P          = atm.mm*atm.mol*const.g/ea.sa;% Pa; sea surface pressure (Pa = kg/m*s2) = 1 bar
atm.N2pal      = atm.mol.*0.78;    % mol; atm n2 present level
atm.H2pal      = 5.5e-5*atm.mol;   % mol; atm H2 present atmospheric level (PAL)
atm.O2pal      = 0.21*atm.mol;     % mol; atm O2 present atmospheric level (PAL)
atm.CH4pal     = 700e-9*atm.mol;   % mol; atm CH4 pre-industrial level (1803 ppb in 2011 is 150% larger than preindustrial est, IPCC Ar 5)
atm.CO2pal     = 300e-6*atm.mol;   % mol; atm CO2 pre-industrial level (280-400 ppm) (Siegenthaler + Sarmiento, 1993)
atm.NH3pal     = 6.7e-10.*atm.mol; % mol; atm NH3 level (Byrne + Goldblatt, 2014, table 1)
atm.PAN        = atm.N2pal.*2;     % mol; atm present level of N
atm.fCO2       = 347.9e-6;         % atm; CO2 fugacity (Zeebe & Wolf-Gladrow 2001)
atm.z          = 12e3;             % m; approximate height of tropopause
atm.mass       = 5.1e18;           % kg; mass of entire atmosphere
atm.trop       = 0.9*atm.mass;     % kg; mass of troposphere (90% of atm)
atm.vol        = atm.z*(ea.rad^2); % m3; volume of troposphere
atm.L          = atm.vol*conv.L2m3;% L; volume of troposphere
atm.rho        = atm.mass/atm.vol; % kg/m3; density of troposphere
atm.O3.lim     = 3e-3;             % PAL; threshold for permanent ozone layer and UV shielding (Catling + Zahnle, 2020; Zahnle et al. 2006)  
atm.K0.N2      = 6.4e-6;           % mol/m3Pa; reference Kh at T0 for N2 (Sander, 2015)
atm.vH.N2      = 1300;             % K; temperature parameter for N2 (Sander, 2015)
atm.K0.NH3     = 5.9e-1;           % mol/m3Pa; reference Kh at T0 for NH3 (Sander, 2015)
atm.vH.NH3     = 4200;             % K; temperature parameter for NH3 (Sander, 2015)
atm.K0.CO2     = 3.3e-4;           % mol/m^3Pa; reference Kh at T0 for CO2 (Sander, 2015)
atm.vH.CO2     = 2400;             % K; Henry solubility Temp dependence for CO2 (Sander, 2015)
atm.K0.CH4     = 1.4e-5;           % mol/m^3Pa; reference Kh at T0 for CH4 (Sander, 2015)
atm.vH.CH4     = 1900;             % K; Henry solubility Temp dependence for CH4 (Sander, 2015)
atm.K0.O2      = 1.2e-5;           % mol/m^3Pa; reference Kh at T0 for O2 (Sander, 2015)
atm.vH.O2      = 1700;             % K; van't Hoff parameter for O2 (Sander, 2015)
atm.D.N2       = 2e-5;             % cm2/s; diffusion constant N2 (CRC Handboook, 98th Ed.)
atm.D.NH3      = 1.5e-5;           % cm2/s; diffusion constant NH3 (CRC Handboook, 98th Ed.)
atm.D.CO2      = 1.45e-5;          % cm2/s; diffusion constant CO2 (CRC Handboook, 98th Ed.)
atm.D.CH4      = 1.55e-5;          % cm2/s; diffusion constant CH4 (CRC Handboook, 98th Ed.)
atm.D.O2       = 1.78e-5;          % cm2/s; diffusion constant O2 (CRC Handboook, 98th Ed.)
atm.kesc       = 2.5e13;           % H2 molecules/cm2 s; constant H-escape rate (via Claire et al. 2006, from Walker 1977)
atm.Tx         = 250;              % K; temperature below which water vapor feedback is unnecessary - used for temp param

% ---------------------------- entire earth -------------------------------
ea.csa         = (1/3).*ea.sa;     % m^2; continental surface area
ea.depth       = 40e3;             % m; approximate depth of cont crust today
ea.vol         = ea.csa.*ea.depth; % m^3; continental crust volume    
ea.cp          = 5e5;              % m^2K/Wyr; earth heat capacity
ea.alb         = 0.3;              % earth albedo (~0.3 to 0.22)
ea.tcrustm     = 2.367e22;         % kg; total mass of oceanic + cont crust (Planetary Scientist's Companion)
ea.ccrustm     = 1.522e22;         % kg; mass of continental crust (Planetary Scientist's Companion, assume that 95-100 % is crystalline for simplicity sake)
ea.ocrustm     = 8.450e21;         % kg; mass of ocean crust (Planetary Scientist's Companion, assume that 95-100 % is crystalline for simplicity sake)
ea.crho        = ea.ccrustm./ea.vol;% kg/m^3; approximate density of mostly granite cont crust
ea.Fe          = (0.05./const.mm.Fe).*ea.ccrustm;% mol; total surface Fe ~ continental concentration of Fe oxide (4.09 wt% (Wedepohl 1995) - 7 wt% (McLennan, 2001))
ea.P           = (mean([0.025 0.1]).*1e-2./const.mm.P).*ea.ccrustm;% mol; total surface P ~ continental concentration of phosphorus (~ 0.025 (Ruttenburg, 2003) - 0.1 wt% (Tiessen, 2018; Paytan + McLaughlin 2007; Rudnick + Gao, 2003) - using mean estimate  
ea.C           = 7.5e21;           % mol; total modern surface C ~ total cont org + inorg C, atm/ocean is minor constituent, modern surface is < 1e22 mol total (cont ests from LCR.m and Catling+Claire, 2006; Treatise on Geochem Vol 2 ch 1 ests 100 ppm in mantle ~ 3e22 mol?)
ea.N           = 7e20;             % mol; using ~5 PAN, from estimate total BSE ~ 7+/-4 PAN (Johnson + Goldblatt, 2015)
ea.Si          = 1e22; %5.3e22 % mol SiO3; continental reservoir of weatherable Ca silicates (cont crust est 66 wt% SiO2; Rudnick + Gao 2003)

% ---------------------------- upper mantle -------------------------------
m.mass         = 1e24;             % kg; upper mantle mass 
m.N            = 2.*atm.PAN;       % mol N; mantle N (mantle res estimated to 2-5x modern atm, Johnson + Goldblatt, 2015)
m.C            = (100e-6.*4e24./const.mm.C);% mol C; 100 ppm in mantle with mass 
m.Fe           = 1.5.*ea.Fe;       % mol Fe; mantle Fe (mantle fluxes approximate a 2 fold change, so we will start with <2x surface total because also outgassing FeO) 
m.Si           = 3.*ea.Si;         % mol SiO3; also approx 2 fold change; tuned!
m.P            = 1.5.*ea.P;        % mol P; start with all earth's P in mantle
% ratios specific to how many pulses we use in the mantle history compensating for the difference in total output vs a purely linear outgassing regime
% for 2 pulses
m.ratio2.C     = 3.578; 
m.ratio2.N     = 3.556; 
m.ratio2.P     = 3.522; 
m.ratio2.SiO3  = 3.513;
m.ratio2.Fe2SiO4= 3.503;
m.ratio2.FeO   = 3.2914; 
% for 3 pulses (depends heavily on timing and the size of mantle reservoir!)
m.ratio3.C     = 2.686; 
m.ratio3.N     = 2.695; 
m.ratio3.P     = 6.0614;
m.ratio3.SiO3  = 2.2041;
m.ratio3.Fe2SiO4= 4.0409;
m.ratio3.FeO   = 2.592; 
m.timescale    = 0:1e6:4.5e9;      % timescale for curve fit
m.pulses       = MantlePulses(3,m);% 3 mantle pulses timed (made in subfunction)
m.Tm0          = 1350+273;         % K; modern upper mantle temp (Korenaga 2010; Padhi et al. 2012; Johnson+Goldblatt 2018)

% ------------------------------- ocean -----------------------------------
oc.S.s         = 35;               % average ocean surface salinity
oc.S.d         = oc.S.s - 1;       % permil; salinity of the deep ocean (WOCE atlas)
oc.S.z         = oc.S.d;           % permil; salinity of deep sediment porewater
oc.S.n         = oc.S.s;           % permil; salinity of shallow sediments
bs = {'s','d','n','z'};
for ib = 1:length(bs)
oc.B.(bs{ib})  = 4.16e-4.*(oc.S.(bs{ib})/35);% mol/kg; total surface ocean boron concentration (Zeebe + Wolf-Gladrow, 2001: A.7.14, Millero et al. 1995)
oc.Ca.(bs{ib}) = 1.028e-2*(oc.S.(bs{ib})/35);% mol/kg; surface ocean [Ca] (Millero 1995, 1982)
end
oc.sa          = (2/3).*ea.sa;     % m^2; surface area of ocean
oc.depth.s     = 100;              % m; depth of surface ocean
oc.depth.d     = 4000;             % m; depth of deep ocean
oc.press.s     = 0.5.*oc.depth.s.*conv.depthbar;% bar pressure of surface ocean, given mid-depth (1 bar = 1e-5 kg/ms2 = 1e-5 N/m2)
oc.press.d     = 0.5.*oc.depth.d.*conv.depthbar;% bar; pressure of deep ocean, given mid-depth 
oc.vol.s       = oc.depth.s*oc.sa; % m^3; surface ocean volume 
oc.vol.d       = oc.depth.d*oc.sa; % m^3; deep ocean volume 
oc.vol.t       = oc.vol.s+oc.vol.d;% m^3; total ocean volume (conservative!)
oc.stag        = 3.96e-5;          % m; stagnant layer boundary thickness (Liss & Slater, 1974)
oc.rho         = 1027;             % kg/m3; density of seawater
oc.mix.sd      = oc.vol.d/1000;    % m^3/yr; volumetric flow rate set by deep ocean reservoir
oc.m.s         = oc.rho*oc.vol.s;  % kg; mass of surface ocean
oc.m.d         = oc.rho*oc.vol.d;  % kg; mass of deep ocean 
oc.m.t         = oc.m.s+oc.m.d;    % kg; mass of entire ocean
oc.hcapJs      = 4200;             % J/kgK; specific heat capacity of water in joules
oc.hcap        = oc.hcapJs/conv.spyr;% Wyr/kgK; specific heat capacity of water in Watts per deg. Kelvin
oc.mixlayer    = oc.sa.*1e5;       % kg; mass of ocean mixed layer
oc.mix.sd      = oc.vol.d./1000;   % m3/yr; volumetric flow rate of ocean boxes, given volume of deep ocean
oc.refPO4      = 0.1e-6.*oc.rho;   % mol/m3; reference concentration of P in ocean when Redfield ratio is 106:1 (0.1 uM, Tanioka + Matsumoto, 2017)
oc.refkrw      = 1.01e-19.*(7^22.4); % /yr; reference rate constant for reverse weathering, given neutral pH
oc.refSi0      = 2.02.^(-5.57.*7); % mol Si/L; reference silica concentration, given neutral pH
oc.refFe       = 100e-6;           % mol/kg; estimated concentrations for Fe(II) in Archean ocean basins (~10-7000 uM, Swanner et al. 2020) 
oc.crustprod   = 3.4e6;            % m2/yr; modern day global average crust production rate (3.4 km2/yr; White + Klein, 2014)
oc.spread      = mean([4,8]).*1e-2;% m/yr; modern day intermediate-rate spreading (4-8 cm/yr; White + Klein, 2014)
oc.dCO2pl      = 4e-5;             % mol/kg; average [CO2] concentration in modern deep ocean/seds
oc.dHCO3pl     = 0.0019;           % mol/kg; average [HCO3] concentration in modern deep ocean/seds
oc.dP.s        = mean([0.1,0.5,1.6]).*1e-6; % mol/kg; average [PO4] in surface ocean (WOCE atlas dataset, pacific/atlantic) 
oc.dP.n        = mean([0.5,1.5,2]).*1e-6; % mol/kg; average [PO4] in shallow seds  (WOCE atlas dataset, pacific/atlantic) 
oc.dP.d        = mean([2,2.5,3]).*1e-6; % mol/kg; average [PO4] in deep ocean (WOCE atlas dataset, pacific/atlantic) 
oc.dP.z        = mean([1.5,2,2.7]).*1e-6; % mol/kg; average [PO4] in pelagic seds (WOCE atlas dataset, pacific/atlantic) 
oc.dN.n        = mean([0.01 0.12]).*1e-5;% mol/kg; estimated range of [NH4] in surface ocean/shelf seds (Gruber, 2008)
oc.dN.z        = mean([0.1 0.55]).*1e-6;% mol/kg; estimated range of [NH4] in deep ocean/seds (Gruber, 2008)
oc.dSi         = 0.1e-3.*conv.L2m3;% mol/m3; upper estimate modern marine dissolved Silica concentration (<= 0.1 mM; Isson + Planavsky, 2018)

% ----------------------------- sediments ---------------------------------
% ** NOTE : all diffusion coefficients assume t = 25C ** 
sed.rho        = 2700;             % kg/m3; sediment density (Johnson + Goldblatt, 2017)
sed.sa.n       = oc.sa.*0.08;      % m^2; surface area of continental shelves (5-10%; Yool + Fasham, 2001)
sed.sa.z       = oc.sa.*0.92;      % m^2; surface area of pelagic sediments, ie. the rest of the ocean
sed.depth.u    = 500;              % m; unreactive sediment depth (Johnson + Goldblatt, 2017)
sed.depth.n    = 0.1;              % m; reactive sediment depth (Johnson + Goldblatt, 2017)
sed.depth.z    = 0.1;              % m; reactive sediment depth (Johnson + Goldblatt, 2017)
oc.press.z     = 2.*oc.press.d;    % bar; pressure of pelagic sediments = full-depth deep ocean
oc.press.n     = 2.*oc.press.s;    % bar; pressure of neritic sediments = full-depth surface ocean
sed.vol.u      = oc.sa*sed.depth.u;% m3; unreactive sediment volume
sed.pore       = 0.7;              % sediment porosity (intermediate between LOSCAR estimates for pure clay ~0.85 and pure carbonate ~0.62)
bnz            = {'n','z'}; 
for ix = 1:length(bnz)
sed.depth.(bnz{ix}) = 0.1;         % m; reactive sediment depth (Johnson + Goldblatt, 2017)
sed.m.(bnz{ix})     = sed.rho*(sed.sa.(bnz{ix})*sed.depth.(bnz{ix}));% kg; sediment mass
sed.vol.(bnz{ix})   = sed.sa.(bnz{ix})*sed.depth.(bnz{ix}); % m3; sediments volume
% porespace water in sediments
oc.depth.(bnz{ix})  = sed.depth.(bnz{ix});        % m; depth of sed porewater = depth of sed
oc.vol.(bnz{ix})    = sed.vol.(bnz{ix}).*sed.pore;% m3; sed porewater volume
oc.m.(bnz{ix})      = oc.rho*oc.vol.(bnz{ix});    % kg; mass of reactive sed layer water 

end
sed.difflength = 0.01;             % m ; characteristic diffusion length (aka L)
sed.diff.HNO3  = 5e-6.*conv.diff;  % m2/yr; NO3 diffusion rate (Butcher et al. 1992 Global Biogeochem Cycles)
sed.diff.NH4   = 7e-6*conv.diff;   % m2/yr; NH4 diffusion rate (Butcher et al. 1992 Global Biogeochem Cycles)
sed.diff.H3PO4 = 3e-6*conv.diff;   % m2/yr; PO4 diffusion rate (Butcher et al. 1992 Global Biogeochem Cycles)
sed.diff.O2    = 2.42e-5*conv.diff;% m2/yr; O2 diffusion rate (CRC handbook online)
sed.diff.NH3   = 1.5e-5*conv.diff; % m2/yr; NH3 diffusion rate (CRC handbook online, Boudreau 1996 = 2.28e-5)
sed.diff.CO2   = 1.91e-5*conv.diff;% m2/yr; CO2 diffusion rate (CRC handbook online, CRC 2017 = 1.45e-5)
sed.diff.N2    = 2e-5*conv.diff;   % m2/yr; N2 diffusion rate (CRC handbook online, Boudreau 1996)  
sed.diff.HCO3  = 5.1e-6*conv.diff; % m2/yr; HCO3 diffusion rate (Bauer et al. 1995, Zeebe 2011 = 1.19e-9 m2/s)
sed.diff.CH4   = 1.84e-5*conv.diff;% m2/yr; CH4 diffusion rate (CRC handbook online, Iversen + Jorgensen, 1993 = 0.87e-5)
sed.diff.CO3   = 0.92e-9*conv.spyr;% m2/yr; CO3 diffusion rate (Zeebe 2011 = 0.92e-9 m2/s)
sed.diff.OH    = 0.5.*sed.diff.O2; % m2/yr; * NO SOURCE, PLACEHOLDER! * just assuming OH is more efficient at diffusing than O2
sed.diff.H     = 5.11e-5*conv.diff;% m2/yr; H2 diffusion rate (CRC online handbook)
sed.diff.RN    = mean([sed.diff.NH3,...
    sed.diff.NH4]);                % m2/yr; reduced N diffusion = mean of NH4/NH3 rates
sed.diff.DIC   = mean([sed.diff.CO2,...
    sed.diff.HCO3,sed.diff.CO3]);  % m2/yr; dissolved inorganic carbon diffusion
sed.diff.TA    = sed.diff.HCO3;    % m2/yr; alkalinity diffusion

% -------------------------- solar constants ------------------------------
S.Pref         = 1361;             % W/m2; present day solar constant (Fs)
S.Aref         = (1-(0.38.*-4e9./4.55e9)).^(-1).*S.Pref;% W/m2; solar constant during the Eoarchean, 4 Ga (using Gough 1981 formulation)
S.Href         = (1-(0.38.*-4.5e9./4.55e9)).^(-1).*S.Pref;% W/m2; solar constant during the Hadean, 4.5 Ga (using Gough 1981 formulation)

% ----------------------------- fractions ---------------------------------
f.Ashelf       = sed.sa.n./oc.sa;  % fraction of shelf surface area to total ocean surface area 
f.Pshelf       = 0.3;	           % most of primary production occurs on continental shelves (20-50%; Longhurst et al 1995; Yool & Fasham 2001)
f.acc          = 0.05;             % fraction of slab accreted (ophiolite emplacement) off at arc setting (seds = f_acc x2)
f.melt         = 0.85;             % fraction of slab/seds melted at arc setting (either volatilized or recrystallized!)
f.sub          = 0.05;             % fraction of sediments subducted at arc setting (slab = f_sub x 2)
f.UMgas        = 0.1;              % assume that 10% of upper mantle is outgassable at any given time       
f.UMred        = 0.01;             % assume that only 1% of outgassed carbon/nitrogen is reduced
% species in slab have different affinities for volatilization/recrystallization (f_cryst = f_melt - f_volc)
sps = {'CaCO3','NH4','SP','CP','FeOH3','Fe2SiO4','ON','OP','OC'}; 
frs = [1, 0.25, 0, 0, 0, 0, 0.25, 0, 1]; % fraction of melt that is volatilized
for ifp = 1:length(sps)
    f.volc.(sps{ifp}) = f.melt .* frs(ifp); 
end
f.scavP        = 0.1;              % maximum scavengable fraction of OP - tuned
f.apatite      = 0.01;             % tuned to get ~ modern inorg P burial ~ 1-5e11; fraction of P (apatite) in sedimentary carbonates, phosphorites (Knudsen + Gunter, 2002- ordovician limestone 5-6 wt% apatite)

% -------------------- half-saturation sensitivities ----------------------
% Half-saturation concentrations for uptake (K_sens)
Ksens.On = 20e-6;                  % mol/L; O2 sensivity in nitr (Fennel et al. 2005) - half sat 
Ksens.Nn = 5.14e-6./const.mm.NH4;  % mol/L; NH4 sensivity in nitrification (from 5.14 mg/L, Dincer + Kargi, 2000)
Ksens.Nd = 0.1e-6./const.mm.NO3;   % mol/L; NO3 sensivity in denitrification (from 0.1-0.3 mg/L, Dincer + Kargi, 2000) -  use lower threshold
Ksens.Oor= 12e-6;                  % mol/L; O2 sensivity in remineralization (from 4-12e-6 mol/L, Laufkotter et al. 2017) - use higher threshold
Ksens.Omt = 5.7e-6;                % mol/L; O2 sensivity methanotrophy (from 5.7 ÂµM, Ren et al. 1997)
Ksens.Hrw = 10^(-7);               % mol/L; sensitivity [H] in reverse weathering (ie. at pH = 7) 
Ksens.FeII= 400e-6;                % mol/L; assumed shelf/surface ocean [Fe(II)] in ferruginous Archean ocean (> 20 - 928 ÂµM; Swann et al. 2020, table 1)

% Half-saturation concentrations for inhibition (K_inhib)
Kin.Nm = 10e-6;                    % mol/L; NO3 inhibiting methanogenesis (anoxic remin; Van Cappellen + Wang, 1995)
Kin.Od = 205e-9;                   % mol/L; O2 inhibiting denitrification (Tiano + Dalsgaard et al, 2014)
Kin.Oar= 8e-6;                     % mol/L; O2 inhibiting anoxic remin (Van Cappellen and Wang, 1995) -- functionally equivalent to Ksens.Oor
Kin.Hn = 10^(-8)  ;                % mol/L; inhibition [H] in nitrification (8Â±0.5 pH, Dincer + Kargi, 2000)
Kin.Hd = 10^(-7.5);                % mol/L; H inhibiting denitrification (from optimal denit 7-8 pH, Dincer + Kargi, 2000)
Kin.Oap= 0.2e-6./const.mm.O2;      % mol/L; O2 inhibiting photoferrotrophy (assuming that this can only occur in anoxic settings, [O2] < 0.2 mg/L; Crowe et al. 2008)

% convert to mol/m3 
ks = fieldnames(Ksens); ki = fieldnames(Kin); 
for is = 1:length(ks)
    Ks.(ks{is}) = Ksens.(ks{is}) .* conv.L2m3; % mol/m3; uptake concentration
end
for in = 1:length(ki)
    Ki.(ki{in}) = Kin.(ki{in}) .* conv.L2m3; % mol/m3; inhibiting concentration
end

% distribution coefficient ([X]_absorbed / [X]_solution)
KdP = 0.07e-6 .* conv.L2m3;        % m3/mol P; distribution coefficient for Feoxide P sorption (0.07 uM-1 == 0.07e-6 kg/mol; Bjerrum + Canfield 2002- backcalculated from BIFs)

% ----------------------------- timescales --------------------------------
taunit         = (1/1.15) ./ 365 ; % yr; nitrification timescale (from k = 1.15 /day ; Dincer + Kargi, 2000)
tauoxrm        = (1/0.0302)./365 ; % yr; oxic remin rate (r = 0.0302 /day ; Kriest + Oschlies, 2008)

% ------------------------------ speciation -------------------------------
% Coefficients for pressure effect on dissociation constants 
% .HM = (From Roberta Hamme's "carbonate_eq7.m" code, based on Millero 1995/1979, and etc. internal references)
spx.A0.HM.k1 = -25.5;              % K1 = H2CO3 dissociation constants (Millero 1979)
spx.A1.HM.k1 = -0.151;                 
spx.A2.HM.k1 = 0.1271; 
spx.A3.HM.k1 = 0; 
spx.B0.HM.k1 = -0.00308;
spx.B1.HM.k1 = -5.78e-4;  
spx.B2.HM.k1 = 0.0000877;        
spx.A0.HM.k2 = -15.82;             % K2 = HCO3 dissociation constants (Millero 1979)
spx.A1.HM.k2 = 0.321;  
spx.A2.HM.k2 = -0.0219;    
spx.A3.HM.k2 = 0; 
spx.B0.HM.k2 = 0.00113;      
spx.B1.HM.k2 = -3.14e-4;
spx.B2.HM.k2 = -0.0001475; 
spx.A0.HM.kb = -29.48;             % Kb = Boric acid dissociation constants (Millero 1979)
spx.A1.HM.kb = 0.295; 
spx.A2.HM.kb = 0.1622;
spx.A3.HM.kb = -0.002608;  
spx.B0.HM.kb = -0.00284;      
spx.B1.HM.kb = 3.54e-4;
spx.B2.HM.kb = 0; 
spx.A0.HM.kw = -20.02;             % Kw = water dissociation constants (Millero 1983)
spx.A1.HM.kw = 0.1119;  
spx.A2.HM.kw = -0.001409;  
spx.B0.HM.kw = -5.13e-3;     
spx.B1.HM.kw = 0.0794e-3;         
spx.A0.HM.kc = -48.76;             % Kc = Calcite dissociation constants (Millero 1983)
spx.A1.HM.kc = 0.5304; 
spx.A2.HM.kc = 0;       
spx.B0.HM.kc = -11.76e-3;    
spx.B1.HM.kc = 0.3692e-3;    
spx.A0.HM.ka = -46;                % Ka = Aragonite dissociation constants (Millero 1983)
spx.A1.HM.ka = 0.5304; 
spx.A2.HM.ka = 0;       
spx.B0.HM.ka = -11.76e-3;   
spx.B1.HM.ka = 0.3692e-3;   
% .ZWG = (From Zeebe + Wolf-Gladrow (2001) and their accompanying "equic.m" code, based on Millero 1995, and etc. internal references)
% code includes errata re: table A.11.1 --> a2/b0/b1 need to be multiplied by 1e-3, not 10e3! And note that coeff. signs disagree in Millero 95/79 
spx.A0.ZWG.k1 = -25.5;             % K1 = H2CO3 dissociation constants 
spx.A1.ZWG.k1 = 0.1271;                 
spx.A2.ZWG.k1 = 0.0e-3; 
spx.B0.ZWG.k1 = -3.08e-3;
spx.B1.ZWG.k1 = 0.0877e-3;  
spx.B2.ZWG.k1 = 0.0;  
spx.A0.ZWG.k2 = -15.82;            % K2 = HCO3 dissociation constants 
spx.A1.ZWG.k2= -0.0219;  
spx.A2.ZWG.k2 = 0.0e-3;    
spx.B0.ZWG.k2 = 1.13e-3;      
spx.B1.ZWG.k2 = -0.1475e-3;
spx.B2.ZWG.k2 = 0.0; 
spx.A0.ZWG.kb = -29.48;            % Kb = Boric acid dissociation constants 
spx.A1.ZWG.kb = 0.1622; 
spx.A2.ZWG.kb = -2.608e-3;
spx.B0.ZWG.kb = -2.84e-3;      
spx.B1.ZWG.kb = 0.0e-3;
spx.B2.ZWG.kb = 0.0 ; 
spx.A0.ZWG.kw = -20.02;            % Kw = water dissociation constants 
spx.A1.ZWG.kw = 0.1119;  
spx.A2.ZWG.kw = -1.409e-3;  
spx.B0.ZWG.kw = -5.13e-3;     
spx.B1.ZWG.kw = 0.0794e-3; 
spx.B2.ZWG.kw = 0.0; 
spx.A0.ZWG.kc = -48.76;            % Kc = Calcite dissociation constants 
spx.A1.ZWG.kc = 0.5304; 
spx.A2.ZWG.kc = 0.0e-3;       
spx.B0.ZWG.kc = -11.76e-3;    
spx.B1.ZWG.kc = 0.3692e-3;   
spx.B2.ZWG.kc = 0.0; 
spx.A0.ZWG.ka = -46;               % Ka = Aragonite dissociation constants 
spx.A1.ZWG.ka = 0.5304; 
spx.A2.ZWG.ka = 0.0e-3;       
spx.B0.ZWG.ka = -11.76e-3;   
spx.B1.ZWG.ka = 0.3692e-3;   
spx.B2.ZWG.ka = 0.0; 

% --------------------- forcings' time dependencies -----------------------
% define transition timescales for the time-dependent evolutions of
% photosynthesizers, fungi, and land plants
td.transitionphoto = 1e8;           % yrs; photosynthesis
td.transitionfungi = 5e7;           % yrs; fungal organsims 
td.transitionplant = 5e7;           % yrs; vascular land plants (and body size)
% define when the transitions should start (starting model at 4 Ga == 0 yr)
td.initphoto = 4e9 - 3.6e9;         % photosynthesis pop growth starts at 3.6 Ga
td.initfungi = 4e9 - 8e8;           % land colonization by fungi circa 600-700 Ma, after evolving between 820-1200 Ma (Heckman et al. 2001, Lucking et al. 2009) -- enhances P weathering through selective mining
td.initbody  = 4e9 - 6.35e8;        % timed body size growth at Ediacaran (635 Ma)
td.initplant = 4e9 - 4e8;           % land organism pop growth after vascular plant colonization ~ 400 Ma
% calculate slopes given a certain length of transition (ie. slow to fast, 1e7 years - 1e8 years) and the desired start/end points (endpoints = 1)
% using increasing sigmoid == y = 1/be^-kt and decreasing sigmoid == y = be^-kt
td.photok = -log(1/1e12).*(1/td.transitionphoto);  % start at 1e-12, end at 1 
td.fungik = -log(1/10) .* (1/td.transitionfungi);  % start at 1e-12, end at 1 
td.bodyk  = -log(1/5) .* (1/td.transitionplant);   % start at 1/5, end at 1
td.plantk = -log(1/1e12) .* (1/td.transitionplant);% start at 1e-12, end at 1
td.terrk  = -log(1/1e12) .* (1/td.transitionfungi);% start at 1e-12, end at 1


%% ------------------------ Radiative forcings ----------------------------
% Calculate eference RFs for a range of pCO2, pNH3 values
gastype = {'CO2','CH4','N2O','NH3'};
% read radiative forcing data from Byrne & Goldblatt, 2014 (Climates of the Past- Copernicus)
addpath ./ByrneSI/
for ii = 1:length(gastype)
    gasname = gastype{ii}; 
    addpath(['./ByrneSI/',gasname]);                                       % add each folder in Byrne SI to the path
    % inputs names of folders starting with 'conc_' (eg. 'conc_2.2')
    c_file = dir(['ByrneSI/',gastype{ii},'/N2_1bar/']);                    % only looking at 1 bar atm
    c_fileref = {c_file.name};

    for jj = 4:length(c_fileref)                                           % start at 4th entry to ignore '.', '..', and './DS_Store'
        cfile = c_fileref{jj};
        filename = ['ByrneSI/',gasname,'/N2_1bar/',cfile,'/cloud_0_0_0.hrt']; % open the datasets with no cloud feedback
        fid = fopen(filename,'r');
        % check file can be opened
        if fid == -1 
            disp(['File ',filename,' can not be opened.']);
        else 
            % take concentration from filename (given as the neg log of the
            % labelled number, ie. conc_2 = 10^-2)
            num = str2double(strrep(cfile,'conc_',''));
            gs.(gasname).neglogconc(jj-3) = num;                           % concentration as negative log_10

            % read rf data from file
            A = textscan(fid,'%f%f%f%f%f%f%f','headerlines',1); 
            barp  = [A{1}];                                                % first column contains atm pressure in bar 
            Flwdn = [A{3}];                                                % third column contains downward LW radiation
            Flwup = [A{4}];                                                % forth column contains upward LW radiation
            % the file is set up as an atmospheric profile, with increasing
            % pressure going down the list. We calculate radiative forcing
            % as "net incoming radiation at TOA" ie. downward radiation -
            % upward radiation at the approximate top of atmosphere 
            gs.(gasname).rf(jj-3) = Flwdn(1) - Flwup(1);                   % net downward LW radiation at the assumed TOA 
        end
        fclose(fid);
    end
    % convert concentrations
    val = -gs.(gasname).neglogconc; 
    gs.(gasname).conc = 10.^val;                                           % officially a concentration, hopefully with all sig figs
end
% reorder output into column vectors and normalize
[bb] = UnpackBB(gs,const,S,ea,atm); 
const.aRF  = bb.a;                                                           % radiative forcing constants NOT the same as temperature parameterization constants
const.bRF  = bb.b; 

%% generate a universal structure to spread the good word of constants!
v.const    = const; 
v.atm      = atm; 
v.oc       = oc; 
v.conv     = conv; 
v.ea       = ea; 
v.m        = m; 
v.sed      = sed;  
v.bb       = bb; 
v.color    = SpeciesColors(); 
v.S        = S; 
v.f        = f; 
v.Ks       = Ks;
v.Ki       = Ki;  
v.KdP      = KdP; 
v.taunit   = taunit; 
v.tauoxrm  = tauoxrm; 
v.spx      = spx; 
v.td       = td; 

%% -------------------- Temperature Parameterization ----------------------
% Use a function to determine the constants aT, bT, q for temperature
% parameterization 
[v.tp]     = TempParam(v,bb.method); 

%% Create a folder for storing model output if ones does not already exist
v.runfolder='OUTPUT'; % feel free to change the names, I like all caps :)
v.figfolder='FIGURES';
if isfolder(v.runfolder) == 0  
    mkdir(pwd,v.runfolder); 
end
if isfolder(v.figfolder) == 0  
    mkdir(pwd,v.figfolder); 
end
addpath(v.runfolder); 
addpath(v.figfolder);

%% subfunction for unpacking and normalizing the Byrne gas dataset
function bb = UnpackBB(gs,const,S,ea,atm)
    % reorder output into column vectors
    gas = {'CO2','CH4','NH3'}; 
    for ig = 1:length(gas)
        rfx = gs.(gas{ig}).rf';                     % W/m2; radiative forcing
        bb.rf.(gas{ig}) = rfx-rfx(end);             % normalize the dataset!
        bb.c.(gas{ig}) = gs.(gas{ig}).conc';        % n/natm concentration        
    end

    %% uncomment to plot the relationships between gas concentrations and RF
%     Plot_GasRFRelationship(bb,atm); 

    %% solve aRF^b for a and b! - NOTE THESE ARE NOT THE SAME as aT and bT for the temperature parameterization
    % these constants equilibrate the radiative forcings for different
    % solar constants, not the effect of the greenhouse gases on
    % temperature!!
    
    % dT/dt = solarin - (TIR out - GHE) where solarin = S/4(1-alb)
    % and TIRout = bol*T^4, GHE = aRF^b
    % aRF^b = bol*T^4 - S/4*(1-alb) 
    era = {'mod','arc'};
    for ie = 1:length(era)
        switch era{ie}
            case 'mod'
                sun = S.Pref;
            otherwise
                sun = S.Aref; % solar constant in Eoarchean
        end
        dif.(era{ie}) = (const.bol.*(289^4) - sun./4.*(1-ea.alb));  % difference between in-out radiation 
    end
    co2pal = atm.CO2pal./atm.mol;
    % METHOD A:
    % first look up the RFs corresponding to 300 ppm and x300ppm for archean
    palCO2.A = 150; % reading off of figure 4a (Goldblatt, McDonald, and McCusker, 2021)
    rfm.A = RFInterp(co2pal,'CO2',bb);
    rfa.A = RFInterp(palCO2.A*co2pal,'CO2',bb);
    
    % METHOD B:
    % assume that the minimum difference between the archean and modern
    % earth temperatures is 30-50 W/m2, such that RF_mod = 35 W/m2 @ 280 ppm
    % CO2, and that RF_arc = RF_mod + (30-50) W/m2. Still use Byrne's data to
    % find the appropriate RF for the modern level of CO2. 
    rfm.B = rfm.A; % ~ 33 W/m^2 at 280 ppm
    rfa.B = rfm.B + 30; 
    Bco2a = find(bb.rf.CO2>(rfa.B-5) & bb.rf.CO2<(rfa.B+5)); %find the [CO2] corresponding to this RF
    palCO2.B = bb.c.CO2(Bco2a(end)).*atm.mol./atm.CO2pal; 
    
    %% Use RFs to find a and b forcing constants, choosing either method A or B
    % a = arc_dif/(rfa^b); 
    % arc_dif/(rfa^b)(rfm^b) = mod_dif --> mod_dif/arc_dif = rfm^b/rfa^b -->
    % mod_dif/arc_dif = (rfm/rfa)^b --> b = ln(mod_dif/arc_dif)/ln(rfm/rfa)
    bb.method = 'A';%'B';% 
    bb.aPAL = palCO2.(bb.method);                                        % archean PAL CO2
    bb.b = log(dif.mod./dif.arc)./log(rfm.(bb.method)./rfa.(bb.method)); % constant b in aRF^b
    bb.a = dif.arc./(rfa.(bb.method).^bb.b);                             % constant a in aRF^b

end

%% Subfunction: plot RF and gas concentrations' relationship
function Plot_GasRFRelationship(bb,atm)
    % quick plot of the [GHG] and RF relationship
    figure(); clf; hold on; 
    gases = {'CO2','CH4','NH3'}; 
    for ig = 1:length(gases)
        gasdisp{ig} = [gases{ig}(1:end-1),'_',gases{ig}(end)]; 
        plot(bb.c.(gases{ig}),bb.rf.(gases{ig}),'DisplayName',gasdisp{ig}); 
    end
    xl = xline(1e-6,'--','DisplayName','1 ppmv'); xl.LineWidth = 2.5; 
    legend('-DynamicLegend','location','northwest'); 
    xlabel('Gas abundance (mol_{gas} / mol_{atm})','FontSize',14); 
    ylabel('Greenhouse Forcing (W/m^2)','FontSize',14); 
    set(gca,'xscale','log'); box on; grid on; 
    
    % and show how these profiles fit with modern GHG levels and RFs
    rfrefs = [33.5,1.4,0.7]; cols = {'b','r','y'}; 
    ccrefs = [atm.CO2pal,atm.CH4pal,atm.NH3pal]./atm.mol;
    for ir = 1:length(gases)
        plot(ccrefs(ir),rfrefs(ir),['o',cols{ir}],'DisplayName',[gasdisp{ir},' modern ref']);
    end
    set(gca,'ylim',[0 40]); 
    
    % and in log!
    figure(); clf; hold on; 
    for is = 1:length(gases)
        plot(-log10(bb.c.(gases{is})),bb.rf.(gases{is}),'DisplayName',gasdisp{is});
    end
    legend('-DynamicLegend','location','northwest'); 
    xlabel('Log Gas abundance -log_{10}(mol_{gas} / mol_{atm})','FontSize',14); 
    ylabel('Greenhouse Forcing (W/m^2)','FontSize',14); 
    set(gca,'xscale','lin','ylim',[0 55],'xlim',[1 9]); box on; grid on; 

    % a plot to highlight the [GHG] that gets us at least 1 W/m2
    figure(); clf; hold on; 
    for iv = 1:length(gases)
        plot(bb.c.(gases{iv}),bb.rf.(gases{iv}),'DisplayName',gasdisp{iv}); 
    end
    yl = yline(1,'--','DisplayName','1 W/m^2'); yl.LineWidth = 2.5; 
    legend('-DynamicLegend','location','northwest'); 
    xlabel('Gas abundance (mol_{gas} / mol_{atm})','FontSize',14); 
    ylabel('Greenhouse Forcing (W/m^2)','FontSize',14); 
    set(gca,'xscale','log','ylim',[0 2]); 
    grid on; box on;     

end

%% Subfunction : Assign colors to model species
% optional input to plot how the colors look

function color = SpeciesColors(PLOT)
if ~exist('PLOT','var')
    PLOT = 'off';
end

% Make uniform colors for different species, determined by major element
spc = {'C','DIC','CO2','OC','CaCO3','CH4','TA','LB','HCO3','CO3'}; % greens
spn = {'N','NH3','NH4','HNO3','N2','RN','ON'};                     % reds/oranges
spp = {'P','H3PO4','FePO4','CP','SP','OP'};                        % blues
spf = {'Fe','Fe2SiO4','Fe2O3','FeO','FeOH3'};                      % purples
sph = {'H','H2O','OH','pH'};                                       % teals
spo = {'O','O2'};                                                  % yellow

% make organic species cooler (bluer), inorganic species warmer (more red), 
% and oxidized species more yellow 
% some species use the same color (ie. OC = LB = CH4) 
for ic = 1:length(spc)
   switch spc{ic}
       case {'C','DIC'}
           color.(spc{ic}) = [0 0.5 0]; 
       case {'OC','CH4','LB'}
           color.(spc{ic}) = color.C + [0.05 0.35 0.4];
       case {'CO2','CaCO3','TA'}
           color.(spc{ic}) = color.C + [0.1.*ic 0.05.*ic 0.1]; 
       case 'CO3'
           color.(spc{ic}) = color.CaCO3; 
       case 'HCO3'
           color.(spc{ic}) = color.TA;
   end
end

for in = 1:length(spn)
    switch spn{in}
        case {'N','N2'}
            color.(spn{in}) = [0.85 0.25 0.25];
        case {'ON','NH3'}
            color.(spn{in}) = color.N + [0.05 0.35 0.34]; 
        case {'RN','NH4'}
            color.(spn{in}) = color.N + [0.15 0 0.15];
        case 'HNO3'
            color.(spn{in}) = color.N + [0 0.25 0]; 
    end
end

for ip = 1:length(spp)
    switch spp{ip}
        case {'P','CP'}
            color.(spp{ip}) = [0 0.5 1];
        case {'H3PO4','OP'}
            color.(spp{ip}) = color.P + [0.4 0.45 0]; 
        case {'SP','FePO4'}
            color.(spp{ip}) = color.P + [0.25 0.1 -0.25]; 
    end
end

for ie = 1:length(spf)
    switch spf{ie}
        case {'Fe','FeO'}
            color.(spf{ie}) = [0.5 0.075 0.4];
        case {'Fe2O3','FeOH3'}
            color.(spf{ie}) = color.Fe + [(0.05*ie) 0 0.4]; 
        case 'Fe2SiO4'
            color.(spf{ie}) = color.Fe + [0.15 0 0.1];
    end
end

for ih = 1:length(sph)
   switch sph{ih}
       case {'H','pH'}
           color.(sph{ih}) = [0.15 0.235 1];
       otherwise 
           color.(sph{ih}) = color.H + [0.1.*ih 0 0]; 
   end
end

for io = 1:length(spo)
   color.(spo{io}) = [0.9 0.8 0.125];  
end


% for all other relevant outputs
color.par1     = [0, 0.45, 0.741];       % parameter 1 color
color.par2     = [0.85, 0.325, 0.098];   % parameter 2 color
color.Archean  = [0.7 0.7 0.7];          % Archean model output color
color.Modern   = [0 0 0];                % modern model output color
color.T        = color.par2;             % temperature color
color.fixN     = color.HNO3+[0.15,0.25,0];% fixed N color
color.No       = [1 1 1];

%% CHECK HOW IT LOOKS
switch PLOT
    case 'on'
        figure(); hold on; 
        fs = fieldnames(color);
        xs = 1:10; 
        for iv = 1:length(fs)
            plot(xs,xs+iv,'color',color.(fs{iv}),'displayname',fs{iv}); 
        end
        legend('-dynamiclegend','location','southoutside','numcolumns',7); 

end

end

%% Subfunction : generate timed mantle pulses
function pulses = MantlePulses(xnm,m)
% input is either 2 or 3, for how many pulses to make
    if xnm == 2
        pulse1 = makedist('normal','mu',1e9,'sigma',2.5e8); % gaussian distribution around 1e8 yr with sigma = 2.5e8 == duration 5e8 yr
        pulse2 = makedist('normal','mu',2e9,'sigma',2.5e8); 
        p1 = pdf(pulse1,m.timescale); p2 = pdf(pulse2,m.timescale); 
        pulses = (p1./max(p1)) + (p2./max(p2)) ;            % divide by the maximum value in this synthetic distribution to get a curve from 0-1

    elseif xnm == 3

        pulse1 = makedist('normal','mu',1e8,'sigma',2.5e8); % gaussian distribution around 1e8 yr with sigma = 2.5e8 == duration 5e8 yr
        pulse2 = makedist('normal','mu',2e9,'sigma',2.5e8); 
        pulse3 = makedist('normal','mu',4e9,'sigma',2.5e8); 
        p1 = pdf(pulse1,m.timescale); p2 = pdf(pulse2,m.timescale);  p3 = pdf(pulse3,m.timescale);
        pulses = (p1./max(p1)) + (p2./max(p2)) + (p3./max(p3)); % divide by the maximum value in this synthetic distribution to get a curve from 0-1

    end

end
