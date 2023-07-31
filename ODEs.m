%% ============================= EONS Model =============================== 
% Julia Horne, 2019
%
% ODE equations for each reservoir detailed by summing up the relevant 
% fluxes. Function is called by the ODE15s solver. Outputs the 
% time steps and the reservoir size at each time step. 

% Optional inputs VAR and VARNAME allow user to change a tunable parameter
% in the model. In the basic 'Run_EONS' script, this functionality is used
% for the mantle reductant sensitivity test.

function dy = ODEs(t,y,inp1,v,VAR,VARNAME)
indx = inp1.indx; % pull out the index structure of all field names

%% evaluate the variable names (if they exist) 
if ~exist('VARNAME','var')  &&  ~exist('VAR','var') % if doesn't exist assign a rando value
    VARNAME = 'xxx'; 
    VAR  = 666; 
end

% determine which parameter to switch 
[inp] = DetermineParam(VAR,VARNAME,inp1);

%% compile y input into reservoir structure
for ii = 1:length(y)-1
    r.(indx.res{ii}).(indx.sp{ii}) = y(ii);
end
T = y(end);

%% run flux functions!
[conc]   = Flux_Spec(r,T,inp,v);       
[flux]   = Flux_BGC(t,conc,r,T,inp,v); 
[gasex]  = Flux_AirSea(r,conc,T,v);
[mix]    = Flux_Mix(r,v);
[gf,Fs,~]= Flux_Temp(t,r.a,inp,v); 

%% combine fluxes, sources, sinks, etc. in ordinary differential equations
% ** NOTE: order of dr equations needs to be the same as the r0 values in
% the master code (ie. set in "InitialConditions" script)
% ** ALSO NOTE: Oxygen and dinitrogen reservoirs are tracked in mol O2 and
% mol N2, respectively!
[dout]   = ODE_System(flux,gasex,mix); 

% calculate temperature change
[dT]     = TempODE(T,r.a,gf,Fs,v);

%% recompile output into y structure with same order as y0 and r0
dy = zeros(size(y));
count = 1;
for ii = 1:length(y)-1 
    dy(count) = dout.([indx.sp{ii},'_',indx.res{ii}]); 
    count = count + 1;
end
dy(end) = dT;  

end

%% Subfunction : the actual ODEs

function [d] = ODE_System(flux,gasex,mix) 
%------------------------- atmosphere reservoir ---------------------------

d.NH3_a  = gasex.NH3 + flux.meta.ON + flux.volc.ON + flux.meta.NH4 + flux.volc.NH4...
    + flux.mantle.NH3 - flux.pholys - flux.ammox - flux.Hesc.NH3; 
d.N2_a   = gasex.N2 + 0.5.*(flux.pholys + flux.ammox + flux.Hesc.NH3) + flux.mantle.N2;
d.CO2_a  = gasex.CO2 + flux.mantle.CO2 + flux.meta.CaCO3 + flux.meta.CO2...
    + flux.volc.CO2  + flux.volc.carb + flux.volc.CaCO3 + 0.5.*flux.methox + 0.5.*flux.Hesc.CH4...
    - flux.wthr.carb - flux.wthr.sil - flux.wthr.NH4 - 1.5.*flux.wthr.CP...
    - 3.*flux.wthr.SP - flux.wthr.Feox - (3/8).*flux.pholys - (190/200).*flux.terrprod;
d.O2_a   = gasex.O2 - flux.methox - flux.wthr.oxi - 0.75.*flux.ammox...
    - 0.5.*flux.wthr.Feox + flux.terrprod; 
d.CH4_a  = gasex.CH4 + flux.mantle.CH4 + flux.meta.CH4 + flux.volc.CH4 ...
    + (10/200).*flux.terrprod + (3/8).*flux.pholys - 0.5.*flux.methox - 0.5.*flux.Hesc.CH4 ;
d.H2O_a  = flux.methox + 1.5.*flux.ammox + 0.75.*flux.pholys + 0.5.*(flux.meta.NH4 + flux.volc.NH4)...
    - flux.wthr.carb - 0.5.*(flux.wthr.sil + flux.wthr.NH4) - 1.5.*flux.wthr.CP - 3.*flux.wthr.SP...
    - flux.Hesc.CH4 - flux.wthr.Feox - flux.terrprod;

%------------------------ surface ocean reservoir -------------------------

d.N2_s   = mix.N2.s + flux.denit.N2.s - flux.prod.N2 - gasex.N2 - mix.diff.N2.n;
d.RN_s   = mix.RN.s + flux.wthr.NH4 + flux.wthr.ON + flux.ammon.RN.s + flux.metha.RN.s ...
    + flux.denit.RN.s - flux.prod.RN - flux.ferrotrophy.RN - flux.nitr.s...
    - mix.diff.RN.n - gasex.NH3;
d.HNO3_s = mix.HNO3.s + flux.nitr.s - flux.denit.HNO3.s - flux.prod.HNO3 - mix.diff.HNO3.n;
d.H3PO4_s= mix.H3PO4.s + flux.wthr.PO4 + flux.wthr.OP + flux.ammon.H3PO4.s...
    + flux.denit.H3PO4.s + flux.metha.H3PO4.s - flux.prod.H3PO4 - flux.sorb.s...
    - mix.diff.H3PO4.n;
d.LB_s   = flux.prod.CO2 - flux.death.OC;
d.ON_s   = flux.death.ON - flux.ammon.RN.s - flux.denit.RN.s - flux.metha.RN.s...
    - flux.export.ON - flux.sed.ON.n;
d.OP_s   = flux.death.OP - flux.ammon.H3PO4.s - flux.denit.H3PO4.s...
    - flux.metha.H3PO4.s - flux.export.OP - flux.sed.OP.n ;
d.OC_s   = flux.death.OC - flux.ammon.OC.s - flux.denit.OC.s...
    - flux.metha.OC.s - flux.export.OC - flux.sed.OC.n;
d.CaCO3_s= flux.precip.s - flux.diss.s - flux.export.CaCO3 - flux.sed.CaCO3.n; 
d.TA_s   = mix.TA.s + flux.wthr.sil + 2.*flux.wthr.carb + 3.*flux.wthr.CP...
    + 3.*flux.wthr.SP + 2.*flux.diss.s - 2.*flux.precip.s - mix.diff.TA.n + flux.prod.NH4;
d.DIC_s  = mix.DIC.s + flux.wthr.sil + 2.*flux.wthr.carb + flux.wthr.oxi...
    + flux.wthr.NH4 + 1.5.*flux.wthr.CP + 3.*flux.wthr.SP + flux.wthr.Feox...
    + flux.diss.s + flux.ammon.OC.s + flux.denit.OC.s + flux.metha.CO2.s...
    - flux.prod.CO2 - flux.precip.s - gasex.CO2 + flux.mtrophy.s - mix.diff.DIC.n;
d.CH4_s  = mix.CH4.s + flux.metha.CH4.s - flux.mtrophy.s - mix.diff.CH4.n - gasex.CH4;
d.O2_s   = mix.O2.s + flux.prod.O2 - flux.ammon.O2.s - 2.*flux.nitr.s...
    - gasex.O2 - 2.*flux.mtrophy.s - flux.photox.O2 - mix.diff.O2.n ; 
d.H2O_s  = flux.wthr.oxi + flux.ammon.H2O.s + flux.denit.H2O.s + flux.metha.H2O.s...
    + flux.nitr.s + 2.*flux.mtrophy.s + flux.precip.s + flux.wthr.Feox + 3.*flux.sorb.s...
    - flux.diss.s - flux.prod.H2O - flux.photox.H2O;
d.FeO_s  = mix.FeO.s - flux.photox.FeO - flux.prod.FeO;
d.FeOH3_s= flux.prod.FeO + flux.photox.FeO - flux.export.FeOH3 ... 
    - flux.sed.FeOH3.n - flux.sorb.s;
d.FePO4_s= flux.sorb.s - flux.sed.FePO4.n - flux.export.FePO4; 

%---------------------- neritic sediment reservoir ------------------------

d.N2_n   = mix.diff.N2.n + flux.denit.N2.n;
d.RN_n   = mix.diff.RN.n + flux.ammon.RN.n + flux.metha.RN.n + flux.denit.RN.n...
    - flux.nitr.n; 
d.HNO3_n = mix.diff.HNO3.n + flux.nitr.n - flux.denit.HNO3.n;
d.H3PO4_n= mix.diff.H3PO4.n + flux.ammon.H3PO4.n + flux.denit.H3PO4.n...
    + flux.metha.H3PO4.n + flux.scav.n - flux.sorb.n - flux.burial.PO4.n;
d.ON_n   = flux.sed.ON.n - flux.ammon.RN.n - flux.denit.RN.n - flux.metha.RN.n...
    - flux.burial.ON.n;
d.OP_n   = flux.sed.OP.n - flux.ammon.H3PO4.n - flux.denit.H3PO4.n ...
    - flux.metha.H3PO4.n - flux.scav.n - flux.burial.OP.n;
d.OC_n   = flux.sed.OC.n - flux.ammon.OC.n - flux.denit.OC.n - flux.metha.OC.n...
    - flux.burial.OC.n;
d.CaCO3_n= flux.sed.CaCO3.n + flux.precip.n - flux.diss.n - flux.burial.CaCO3.n...
    - 1.5.*flux.burial.PO4.n;  
d.TA_n   = mix.diff.TA.n + 2.*flux.diss.n - 2.*flux.precip.n...
    - flux.revweather.n - 3.*flux.burial.FePO4.n; 
d.DIC_n  = mix.diff.DIC.n + flux.diss.n + flux.ammon.OC.n + flux.denit.OC.n...
    + flux.metha.CO2.n + flux.mtrophy.n + 1.5.*flux.burial.PO4.n - flux.precip.n; 
d.CH4_n  = mix.diff.CH4.n + flux.metha.CH4.n - flux.mtrophy.n;
d.O2_n   = mix.diff.O2.n - flux.ammon.O2.n - 2.*flux.nitr.n - 2.*flux.mtrophy.n ;
d.H2O_n  = flux.ammon.H2O.n + flux.denit.H2O.n + flux.metha.H2O.n...
    + flux.nitr.n + flux.precip.n - flux.diss.n + 2.*flux.mtrophy.n ...
    + 0.5.*flux.revweather.n + 1.5.*flux.burial.PO4.n + 3.*flux.sorb.n;
d.FeOH3_n= flux.sed.FeOH3.n - flux.sorb.n - flux.burial.FeOH3.n;
d.FePO4_n= flux.sed.FePO4.n + flux.sorb.n - flux.burial.FePO4.n; 

%-------------------------- deep ocean reservoir --------------------------

d.N2_d   = mix.N2.d + flux.denit.N2.d - mix.diff.N2.z;
d.RN_d   = mix.RN.d + flux.ammon.RN.d + flux.metha.RN.d + flux.denit.RN.d...
    - flux.nitr.d - mix.diff.RN.z;
d.HNO3_d = mix.HNO3.d + flux.nitr.d - flux.denit.HNO3.d - mix.diff.HNO3.z;
d.H3PO4_d= mix.H3PO4.d + flux.ammon.H3PO4.d + flux.denit.H3PO4.d + flux.metha.H3PO4.d ...
    - mix.diff.H3PO4.z - flux.sorb.d;
d.ON_d   = flux.export.ON - flux.ammon.RN.d - flux.denit.RN.d - flux.metha.RN.d ...
    - flux.sed.ON.z;
d.OP_d   = flux.export.OP - flux.ammon.H3PO4.d - flux.denit.H3PO4.d - flux.metha.H3PO4.d ...
    - flux.sed.OP.z;
d.OC_d   = flux.export.OC - flux.ammon.OC.d - flux.denit.OC.d - flux.metha.OC.d ...
    - flux.sed.OC.z;
d.CaCO3_d= flux.export.CaCO3 + flux.precip.d - flux.diss.d - flux.sed.CaCO3.z;
d.TA_d   = mix.TA.d + 2.*flux.diss.d - 2.*flux.precip.d - mix.diff.TA.z;
d.DIC_d  = mix.DIC.d + flux.diss.d + flux.ammon.OC.d + flux.denit.OC.d + flux.metha.CO2.d ...
    - flux.precip.d + flux.mtrophy.d - mix.diff.DIC.z;
d.CH4_d  = mix.CH4.d + flux.metha.CH4.d - flux.mtrophy.d - mix.diff.CH4.z;
d.O2_d   = mix.O2.d - flux.ammon.O2.d - 2.*flux.nitr.d - 2.*flux.mtrophy.d...
    - mix.diff.O2.z; 
d.H2O_d  = flux.ammon.H2O.d + flux.denit.H2O.d + flux.metha.H2O.d + flux.nitr.d...
    + 2.*flux.mtrophy.d + flux.precip.d - flux.diss.d + 3.*flux.sorb.d;
d.FeO_d  = mix.FeO.d + flux.mantle.FeO;
d.FeOH3_d= flux.export.FeOH3 - flux.sed.FeOH3.z - flux.sorb.d;
d.FePO4_d= flux.export.FePO4 + flux.sorb.d - flux.sed.FePO4.z; 

%----------------------- reactive sediment reservoir ----------------------

d.N2_z   = mix.diff.N2.z + flux.denit.N2.z; 
d.RN_z   = mix.diff.RN.z + flux.ammon.RN.z + flux.metha.RN.z + flux.denit.RN.z...
    - flux.nitr.z - flux.hyd ;
d.HNO3_z = mix.diff.HNO3.z + flux.nitr.z - flux.denit.HNO3.z;
d.H3PO4_z= mix.diff.H3PO4.z + flux.ammon.H3PO4.z + flux.denit.H3PO4.z ...
    + flux.metha.H3PO4.z + flux.scav.z - flux.sorb.z - flux.burial.PO4.z;
d.ON_z   = flux.sed.ON.z - flux.ammon.RN.z - flux.denit.RN.z - flux.metha.RN.z...
    - flux.burial.ON.z;
d.OP_z   = flux.sed.OP.z - flux.ammon.H3PO4.z - flux.denit.H3PO4.z - flux.metha.H3PO4.z...
    - flux.scav.z - flux.burial.OP.z;
d.OC_z   = flux.sed.OC.z - flux.ammon.OC.z - flux.denit.OC.z - flux.metha.OC.z...
    - flux.burial.OC.z;
d.CaCO3_z= flux.sed.CaCO3.z + flux.precip.z - flux.diss.z - flux.burial.CaCO3.z...
    - 1.5.*flux.burial.PO4.z;
d.TA_z   = mix.diff.TA.z + 2.*flux.diss.z - 2.*flux.precip.z - flux.revweather.z ...
    - 3.*flux.burial.FePO4.z; 
d.DIC_z  = mix.diff.DIC.z + flux.diss.z + flux.ammon.OC.z + flux.denit.OC.z ...
    + flux.metha.CO2.z + flux.mtrophy.z + 1.5.*flux.burial.PO4.z - flux.precip.z - flux.sfw;
d.CH4_z  = mix.diff.CH4.z + flux.metha.CH4.z - flux.mtrophy.z;
d.O2_z   = mix.diff.O2.z - flux.ammon.O2.z - 2.*flux.nitr.z - 2.*flux.mtrophy.z; 
d.H2O_z  = flux.ammon.H2O.z + flux.denit.H2O.z + flux.metha.H2O.z + 0.5.*flux.hyd...
    + flux.nitr.z + flux.precip.z + 2.*flux.mtrophy.z - flux.diss.z ...
    + 0.5.*flux.revweather.z + 1.5.*flux.burial.PO4.z + 3.*flux.sorb.z;
d.FeOH3_z= flux.sed.FeOH3.z - flux.sorb.z - flux.burial.FeOH3.z;
d.FePO4_z= flux.sed.FePO4.z + flux.sorb.z - flux.burial.FePO4.z;

%---------------------- unreactive sediment reservoir ---------------------

d.SP_u   = flux.burial.FePO4.z - flux.acc.SP - flux.cryst.SP - flux.subduct.SP; 
d.CP_u   = flux.burial.PO4.z - flux.acc.CP - flux.cryst.CP - flux.subduct.CP;
d.ON_u   = flux.burial.ON.z - flux.acc.ON - flux.cryst.ON - flux.volc.ON - flux.subduct.ON; 
d.OP_u   = flux.burial.OP.z - flux.acc.OP - flux.cryst.OP - flux.subduct.OP;
d.OC_u   = flux.burial.OC.z - flux.acc.OC - flux.volc.OC - flux.subduct.OC;
d.CaCO3_u= flux.burial.CaCO3.z - flux.acc.carb - flux.volc.carb - flux.subduct.carb; 
d.FeOH3_u= flux.burial.FeOH3.z + flux.burial.FePO4.z - flux.acc.FeOH3...
    - flux.cryst.FeOH3 - flux.subduct.FeOH3;

%----------------------- continental crust reservoir ----------------------

d.NH4_c  = flux.acc.NH4 + flux.cryst.NH4 + flux.cryst.ON - flux.wthr.NH4 - flux.meta.NH4;
d.SP_c   = flux.mantle.P + flux.burial.FePO4.n + flux.acc.SP + flux.cryst.SP ...
    + flux.cryst.OP + flux.cryst.CP - flux.wthr.SP;
d.CP_c   = flux.burial.PO4.n + flux.acc.CP + flux.meta.OP - flux.wthr.CP;
d.ON_c   = flux.burial.ON.n + flux.acc.ON - flux.wthr.ON - flux.meta.ON; 
d.OP_c   = flux.burial.OP.n + flux.acc.OP - flux.wthr.OP - flux.meta.OP; 
d.OC_c   = flux.burial.OC.n + flux.acc.OC + (90/100).*flux.terrprod - flux.wthr.oxi - flux.meta.OC;
d.CaCO3_c= flux.burial.CaCO3.n + flux.acc.CaCO3 + flux.acc.carb - flux.wthr.carb - flux.meta.CaCO3;
d.Fe2SiO4_c= flux.mantle.Fe2SiO4 - flux.wthr.Feox; 
d.FeOH3_c= flux.burial.FeOH3.n + flux.burial.FePO4.n + flux.acc.FeOH3 + flux.cryst.FeOH3; 
d.Fe2O3_c= flux.wthr.Feox; 
d.SiO3_c = flux.mantle.SiO3;

%------------------------- ocean crust reservoir --------------------------

d.NH4_o  = flux.hyd - flux.volc.NH4 - flux.acc.NH4 - flux.cryst.NH4 - flux.subduct.NH4; 
d.CaCO3_o= flux.sfw - flux.volc.CaCO3 - flux.acc.CaCO3 - flux.subduct.CaCO3; 

%------------------------- upper mantle reservoir -------------------------
 
d.N_m    = flux.subduct.NH4 + flux.subduct.ON - flux.mantle.N; 
d.P_m    = flux.subduct.SP + flux.subduct.CP + flux.subduct.OP - flux.mantle.P; 
d.C_m    = flux.subduct.CaCO3 + flux.subduct.carb + flux.subduct.OC - flux.mantle.C;
d.Fe_m   = flux.subduct.FeOH3 - flux.mantle.FeO - 2.*flux.mantle.Fe2SiO4;
d.SiO3_m = - flux.mantle.SiO3;

%--------------------------- Mass Conservation ----------------------------
% add lost untracked species, subtract untracked added species
d.H_x    = flux.Hesc.H + 3.*(flux.meta.OP + flux.cryst.OP)...
    + 4.*flux.subduct.NH4 + 3.*(flux.subduct.ON - flux.mantle.NH3) - flux.cryst.ON... 
    + 2.*flux.subduct.OC - 4.*flux.mantle.CH4...
    + 3.*(flux.subduct.OP + flux.subduct.FeOH3);
d.O_x    = flux.subduct.OC - 2.*flux.mantle.CO2 + 3.*(flux.subduct.CaCO3 + flux.subduct.carb) ... 
    + flux.volc.CaCO3 + flux.volc.carb + flux.meta.CaCO3 ...
    + 0.5.*(flux.revweather.n + flux.revweather.z) - 0.5.*flux.wthr.sil - flux.sfw...
    + 0.5.*(flux.hyd - flux.wthr.NH4 - flux.meta.NH4 - flux.volc.NH4)... 
    + 4.*(flux.subduct.SP + flux.subduct.CP + flux.subduct.OP - flux.mantle.P) ...
    + 3.*flux.subduct.FeOH3 - flux.mantle.FeO - 4.*flux.mantle.Fe2SiO4 + 2.*flux.wthr.Feox; 


end

%% Subfunction : temperature change

function dT = TempODE(T,ra,gf,Fs,v)
Fbb  = v.const.bol.*(T.^4);                % W/m^2; thermal outflux due to blackbody radiation
X    = v.ea.sa./(v.oc.m.t.*v.oc.hcap);     % Km2/Wyr; total ocean heat capacity per year

% calculate pressure broadening (N) effect
N    = (ra.N2./v.atm.N2pal).^v.tp.q; 
% calculate thermal IR opacity (y_ir) with pressure broadening parameter
yir  = N.*(v.tp.aT.*(T-v.atm.Tx) + gf.t.*v.tp.bT);
Fir  = (1 + 0.75.*yir);                    % W/m^2; thermal retention flux from GHG forcing

% Change in temperature from GHG partial pressures + in/outflux radiation 
% balance at current surface 
dT   = X .* (Fs .* Fir  - Fbb); 

end

%% Subfunction : Change parameter value

% You can change a single parameter or multiple at once in a test; enter
% the VAR and VARNAMES as cells if changing multiple. 
function [inpx] = DetermineParam(VAR,NAME,inp)

inpx = inp;             % initialize a new input structure for parameters
if iscell(VAR)          % if the input is a cell, do nothing
    var = VAR;
    name = NAME;
else % if not, assign a dummy cell so that the following switcharoo doesn't send an error message
    var{1} = VAR; 
    name{1} = NAME; 
end
% check the NAME (should match the name of a parameter in 'TunableParameters'
% EXACTLY) and switch to new value VAR
for iv = 1:length(var)
    if ~strcmpi(name{iv},'xxx') % unless this variable is assigned a dummy name, evaluate the rest!
        inpx.(name{iv}) = var{iv}; 

    end
end    
end
