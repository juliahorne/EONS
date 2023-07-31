% Run only the nominal run section of the master script
clear; 
TunableParameters;                          % contains standards and all the reference parameters used in the model 
[inp.indx,y0] = InitialConditions(inp,v);   % initial reservoir sizes, indexed names added to input structure 

% the nominal run assumes photosynthesis arises at 1e9 yr (3.5 Ga) - users
% may choose to input different start/restart times (inp.S1, inp.S2),
% assigned in TunableParameters.m

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


%% Subfunction: Run the EONS model
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
[ot.flux]    = Flux_BGC(tout,ot.conc,ot.r,ot.T,inp,v);             	% function for biogeochemical fluxes 
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

