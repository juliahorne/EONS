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
RunParameterTest(y0,inp,ops,v,tests,values,'fman'); 

% ----------------------------- THE END -----------------------------------

%% Subfunction: Evaluate the model with different mantle influx treatments
function RunParameterTest(y0,inp1,ops,v,tests,varlist,VARNAME)
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
alltests = cell(length(tests),length(varlist));
for ia = 1:length(tests)
    for ib = 1:length(varlist)
        tname = [VARNAME,'_',tests{ia},'_',num2str(varlist(ib))];
        alltests{ia}{ib} = load([v.runfolder,'/',tname,'.mat']);
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
