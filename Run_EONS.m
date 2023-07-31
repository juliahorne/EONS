%% ============================= EONS Model =============================== 
% ----------------------------- Master script -----------------------------
% Julia Horne, 2023

% Check out the README for more detail. We recommend running the following
% sections individually because of long model runtime. 

%% ----------------------- Forward evolution run --------------------------
% The following section produces the model output in the EONS paper (TBD) 
% 'Nominal Run' results.
% We start the model without photosynthesis for 1e9 years, then restart the
% model at that exact point with photosynthesis turned on! Then we splice 
% the runs together.                                    *** footnote 1 ***
    
clear; 
TunableParameters;                          % contains standards and all the reference parameters used in the model 
[inp.indx,y0] = InitialConditions(inp,v);   % initial reservoir sizes, indexed names added to input structure 

% run stage 1 : NO OXYGENIC PHOTOSYNTHESIS 
[t1,y1]       = RunEONSModel(y0,inp.S1,0,inp,ops,v,'S1');  

% run stage 2 : OXYGENIC PHOTOSYNTHESIS
use output as initial conditions for the next stage
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

%% ------------------------ Test Mantle Tuning ----------------------------
% Evaluate how the initial mantle tuning and linear decay effect the GOE
% timing.                                               *** footnote 2 ***
clear; 
TunableParameters;  
[inp.indx,y0] = InitialConditions(inp,v);  
% set up mantle values to test 
values        = 1:1:10; 
tests         = {'constant','linear'}; 
% subfunction executes the tests in parallel
alltests      = RunParameterTest(y0,inp,ops,v,tests,values,'fman');
% Plot the mantle sensitivity tests - many test will take a while to plot
Plot_MantleTest(alltests); 

% ----------------------------- THE END -----------------------------------

%% Subfunction: Run the basic EONS model
function [t,y] = RunEONSModel(y0,STOP,TRESET,inp,ops,v,printname)
    inp.treset = TRESET;                                    % assign reset time 
    [t,y] = ode15s(@ODEs,[0 STOP],y0,ops,inp,v);            % run!
    save([v.runfolder,'/EONS_',printname,'.mat'],'-v7.3');  % save!
end

%% Subfunction: Evaluate the model with different mantle influx treatments
function alltests = RunParameterTest(y0,inp1,ops,v,tests,varlist,VARNAME)
inp1.treset = 0;
STOP1 = inp1.S1; 
STOP2 = inp1.S2;
for it = 1:length(tests)
    TESTNAME = tests{it};                           % either constant or linearly decreasing mantle reductant influx
    inp1.testname = TESTNAME;
    parfor ii = 1:length(varlist)
        TESTNUM = num2str(varlist(ii)); 
        VAR = varlist(ii);                          % the parameter value to test in this run
        inprs = inp1;                               % final input structure
        try
            [t1,y1]       = ode15s(@ODEs,[0 STOP1],y0,ops,inp1,v,VAR,VARNAME);
            inp2          = inp1; 
            inp2.treset   = t1(end); 
            [t2,y2]       = ode15s(@ODEs,[0 STOP2],y1(end,:),ops,inp2,v,VAR,VARNAME);
            % merge the runs for output
            t3            = [t1; t1(end)+t2];
            y3            = [y1; y2]; 
            inprs.(VARNAME)= VAR;                    % remember to reassign or this is overwritten in the output!
            inprs.treset  = 0;                      

        catch ME
            disp(ME.message)
            t3 = zeros(size(y0)); y3 = y0; 
        end
        [output] = CompileOutput(y3,t3,inprs,v);
        parsave(output,v,[VARNAME,'_',TESTNAME,'_',TESTNUM]); 
    end
end
% once all tests are done, save the entire test
clear output
alltests = cell(length(tests));
for ia = 1:length(tests)
    for ib = 1:length(varlist)
        tname = [VARNAME,'_',tests{ia},'_',num2str(varlist(ib))];
        alltests{ia}{ib} =  load([v.runfolder,'/',tname,'.mat']);
        clear output
    end
end
save([v.runfolder,'/',VARNAME,'_Tests.mat'],'alltests','-v7.3');
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

%% Subfunction: save output in parallel run
function parsave(output,vars,testname)
    out = output; 
    v = vars; 
    save([v.runfolder,'/',testname,'.mat'],'-v7.3');
end
 
%% ---------------------------- FOOTNOTES ---------------------------------

%                               *** 1 *** 
% The published (TBD) nominal run assumes photosynthesis arises at 1e9 yr 
% (3.5 Ga) - users may choose to input different start/restart times 
% by changing the values for 'S1' and 'S2' in TunableParameters.m
% Having the 'inp.treset' value as an input to the model is very important;
% otherwise, time sensitive forcings (tdep = the biological transitions and
% the solar constant) will be 'out of time' with the rest of the model run
% after the reset! This input allows the model to change from a 'model
% runtime' = 0 yr at the start of the second run to a 'true time' = 1e9 yr. 
% treset is set in the subfunction. 

%                               *** 2 *** 
% This same setup can be used to test other tunable parameters (assigned in
% the 'TunableParameters.m' script); to do so, one would need to determine 
% a reasonable array of values to test, since the model will likely be very
% sensitive to these parameters and will fail in extreme cases. But that's
% your perogative! The 'tests' input flag only applies to this example, and
% would be extranneous in testing other parameters unless changes are made
% in the flux functions to allow similar linear variations to parameter
% values. If the user wants to try say, linearly decreasing dissolved
% silica concentrations (inp.dSi) then they could use the 'EvolveParam' 
% subfunction in 'Flux_BGC' to do so, but would have to add a statement
% like the one used in assigning the mantle influx. 
