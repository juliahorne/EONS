%% ============================= EONS Model =============================== 
% ---------------------------- Master script -----------------------------
% Julia Horne, 2023

% The following section produces the model output in the EONS paper (TBD) 
% 'Nominal Run' results.
% We start the model without photosynthesis for 1e9 years, then restart the
% model at that exact point with photosynthesis turned on! Then we splice 
% the runs together.                                       *** footnote ***
    
clear; 
TunableParameters;                          % contains standards and all the reference parameters used in the model 
[inp.indx,y0] = InitialConditions(inp,v);   % initial reservoir sizes, indexed names added to input structure 

% run stage 1 : NO OXYGENIC PHOTOSYNTHESIS 
[t1,y1]       = RunEONSModel(y0,inp.S1,0,inp,ops,v,'S1');  

% run stage 2 : OXYGENIC PHOTOSYNTHESIS
% use output as initial conditions for the next stage
[t2,y2]       = RunEONSModel(y1(end,:),inp.S2,t1(end),inp,ops,v,'S2');  

% NOW MERGE the two outputs and compile!
t3            = [t1; t1(end)+t2]; 
y3            = [y1; y2]; 
inp.treset    = 0;
[out]         = CompileOutput(y3,t3,inp,v);
save([v.runfolder,'/','EONS_Full.mat'],'-v7.3');

% Recreate nominal results figures
Plot_NominalRun(out,v);      
Plot_NutrientLimitations(out,v); 
Plot_OxygenHistory(out,v); 

% ----------------------------- THE END -----------------------------------

%% Subfunction: Run the basic EONS model
function [t,y] = RunEONSModel(y0,STOP,TRESET,inp,ops,v,printname)
    inp.treset = TRESET;                                    % assign reset time 
    [t,y] = ode15s(@ODEs,[0 STOP],y0,ops,inp,v);            % run!
    save([v.runfolder,'/EONS_',printname,'.mat'],'-v7.3');  % save!
end

 
%% Subfunction: Compile model output into fluxes/reservoirs/etc. 
% re-run the model component functions with the output y vector to get all
% of the fluxes/reservoirs/concentrations/etc. through time

function [ot] = CompileOutput(yout,tout,inp,v)

% recreate 'r' reservoir structure using indexed reservoirs/species and the
% y output from the ODE solver
for iz = 1:length(inp.indx.res)
    ot.r.(inp.indx.res{iz}).(inp.indx.sp{iz}) = yout(:,iz);
end
ot.T = yout(:,end); 

% run model otput through flux functions to plot flux changes with time
[ot.conc]    = Flux_Spec(ot.r,ot.T,inp,v);                              % function for speciation state 
[ot.flux]    = Flux_BGC(tout,ot.conc,ot.r,ot.T,inp,v);                  % function for biogeochemical fluxes 
[ot.gas]     = Flux_AirSea(ot.r,ot.conc,ot.T,v);                        % function for air-sea gas exchange
[ot.mix]     = Flux_Mix(ot.r,v);                                        % function for ocean mixing fluxes
[ot.rf,~,~]  = Flux_Temp(tout,ot.r.a,inp,v);                            % function for equilibrium temperature and radiative forcing

% run output through functions that show time dependent forcings and timescales
ot.tdep      = Forcings(tout,v);
ot.tau       = Timescales(ot.tdep,inp,v); 

% then put the output cells into a single structure for simplicity
ot.inp       = inp; 
ot.indx      = inp.indx; 
ot.y         = yout; 
ot.t         = tout;

end

%% ---------------------------- FOOTNOTES ---------------------------------

%                                *** *** 
% The published (TBD) nominal run assumes photosynthesis arises at 5e8 yr 
% (3.5 Ga) - users may choose to input different start/restart times 
% by changing the values for 'S1' and 'S2' in TunableParameters.m
% Having the 'inp.treset' value as an input to the model is very important;
% otherwise, time sensitive forcings (tdep = the biological transitions and
% the solar constant) will be 'out of time' with the rest of the model run
% after the reset! This input allows the model to change from a 'model
% runtime' = 0 yr at the start of the second run to a 'true time' = 5e8 yr. 
% treset is set in the subfunction. 
