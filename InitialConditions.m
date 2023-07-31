%% ============================= EONS Model =============================== 
% Julia Horne, 2018

% Initializing reservoirs for the Eoarchean Earth (4 Ga)

function [indx,y0] = InitialConditions(ref,v)
% start with everything == 0
SP.a = {'NH3','N2','CO2','O2','CH4','H2O'};  % atmosphere species
SP.s = {'N2','RN','HNO3','H3PO4','LB','ON','OP','OC','CaCO3','TA','DIC','CH4',...
    'O2','H2O','FeO','FeOH3','FePO4'};       % surface ocean species
SP.n = {'N2','RN','HNO3','H3PO4','ON','OP','OC','CaCO3','TA','DIC','CH4',...
    'O2','H2O','FeOH3','FePO4'};             % neritic sediment species
SP.d = {'N2','RN','HNO3','H3PO4','ON','OP','OC','CaCO3','TA','DIC','CH4',...
    'O2','H2O','FeO','FeOH3','FePO4'};       % deep ocean species
SP.z = {'N2','RN','HNO3','H3PO4','ON','OP','OC','CaCO3','TA','DIC','CH4',...
    'O2','H2O','FeOH3','FePO4'};             % pelagic sediment species
SP.u = {'SP','CP','ON','OP','OC','CaCO3','FeOH3'};% unreactive sediment species
SP.c = {'NH4','SP','CP','ON','OP','OC','CaCO3','Fe2SiO4','FeOH3','Fe2O3',...
    'SiO3'};                                 % continental crust species
SP.o = {'NH4','CaCO3'};                      % ocean crust species
SP.m = {'N','P','C','Fe','SiO3'};            % mantle species
SP.x = {'H','O'};                            % mass conservation

indx.rl = fieldnames(SP);                    % list of all reservoir names
for ia = 1:length(indx.rl)
    RES = indx.rl{ia};
    for ib = 1:length(SP.(RES))
       r0.(RES).(SP.(RES){ib}) = 0;          % structure (r0) with 0 for all species' initial reservoir size
    end
end

%% fill in some reservoirs with Eoarchean estimates

T0          = 288;                          % K; temperature
% atmosphere reservoir, field 'a'
r0.a.NH3    = ref.fNH3.*v.atm.mol;          % mol N; atm ammonia 
r0.a.N2     = v.atm.N2pal;                  % mol N; atm dinitrogen
r0.a.CO2    = ref.fCO2.*v.atm.CO2pal;       % mol C; atm CO2 
r0.a.CH4    = v.atm.CH4pal;                 % mol CH4; atm methane
r0.a.H2O    = 1e18;

% surface ocean
r0.s.N2     = 3e17;                         % mol N; dissolved dinitrogen
r0.s.RN     = 1e14;                         % mol N; dissolved ammonium/ammonia
r0.s.H3PO4  = 1.75e-7*v.oc.m.s;             % mol P; dissolved phosphate (archean est nearly 0.2 uM; Bjerrum + Canfield 2002, Fennel et al. 2005)
r0.s.TA     = 70e-3*v.oc.m.s;               % mol C; ocean alkalinity (2.35e-3 == modern)
r0.s.DIC    = 75e-3*v.oc.m.s;               % mol C; dissolved inorganic carbon (~ 2e-3 == modern) 
r0.s.CH4    = 1e6;                          % mol C; dissolved methane
r0.s.H2O    = 1e19; 
r0.s.FeO    = v.oc.refFe.*v.oc.m.s;         % mol Fe; ferrous iron 

% continental shelf reactive sediments (neritic)
r0.n.N2     = 3e17;                         % mol N; dissolved dinitrogen
r0.n.RN     = 1e9;                          % mol N; dissolved ammonium
r0.n.H3PO4  = 1.75e-7*v.oc.m.n;             % mol P; dissolved phosphate 
r0.n.TA     = 70e-3*v.oc.m.n;               % mol C; ocean alkalinity
r0.n.DIC    = 75e-3*v.oc.m.n;               % mol C; dissolved inorganic carbon
r0.n.CH4    = 1e6;                          % mol C; dissolved methane
r0.n.H2O    = 1e18; 

% deep ocean
r0.d.N2     = 1e18;                         % mol N; dissolved dinitrogen
r0.d.RN     = 1e15;                         % mol N; dissolved ammonium/ammonia
r0.d.H3PO4  = 1.75e-7*v.oc.m.d;             % mol P; dissolved phosphate 
r0.d.TA     = 70e-3*v.oc.m.d;               % mol C; sedimentary alkalinity
r0.d.DIC    = 75e-3*v.oc.m.d;               % dissolved inorganic carbon
r0.d.CH4    = 1e6;                          % mol C; dissolved methane
r0.d.H2O    = 1e20; 
r0.d.FeO    = v.oc.refFe.*v.oc.m.d;         % mol Fe; ferrous iron 

% pelagic reactive sediments 
r0.z.N2     = 1e9;                          % mol N; sedimentary dinitrogen       
r0.z.RN     = 1e9;                          % mol N; sedimentary ammonium/ammonia
r0.z.H3PO4  = 1.75e-7*v.oc.m.z;             % mol P; sedimentary phosphate
r0.z.TA     = 70e-3*v.oc.m.z;               % mol C; sedimentary alkalinity
r0.z.DIC    = 75e-3*v.oc.m.z;               % dissolved inorganic carbon
r0.z.CH4    = 1e6;                          % mol C; dissolved methane
r0.z.H2O    = 1e18; 

% continental crust initially 5% modern emergence, mostly Fe and P bearing
% silicates
r0.c.SP     = 0.05.*v.ea.P;                 % mol P; silicate-bound phosphate in continental crust
r0.c.Fe2SiO4= 0.05.*v.ea.Fe;                % mol Fe; ferrous iron in continental crust
r0.c.SiO3   = 0.05.*v.ea.Si;                % mol Si; crustal silicate reservoir CaSiO3

% start with modern surface reservoirs in upper mantle for these species
r0.m.N      = v.m.N;                        % mol N; mantle N (estimated to 2-5x modern atm, Johnson + Goldblatt, 2015)
r0.m.P      = v.m.P;                        % mol P
r0.m.C      = v.m.C;                        % mol C; 100 ppm 
r0.m.Fe     = v.m.Fe;                       % mol Fe
r0.m.SiO3   = v.m.Si;                       % mol Si

%% compile r0 into a non-structure array y0 (ODE solver can not take in a structure)
count = 1;
for ii = 1:length(indx.rl)
    sl = fieldnames(r0.(indx.rl{ii}));      % only need the first letter of each reservoir string
    for jj = 1:length(sl)
        y0(count) = r0.(indx.rl{ii}).(sl{jj});
        % index the reservoir and species names to create structure in ODE function
        indx.res{count} = indx.rl{ii};      % index for reservoir reference
        indx.sp{count} = sl{jj};            % index for species reference
        count = count + 1;

    end
end
y0(end+1) = T0; 

end
