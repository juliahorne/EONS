%% ============================= EONS Model =============================== 
% Julia Horne, 2019

% Unpack the model output 'out' structure for easier plotting/evaluation 

function [t,r,conc,flux,gasex,mix,RF,T,inp,tdep,tau,indx] = UnpackOutput(out)
    r       = out.r;                            % reservoirs
    conc    = out.conc;                         % concentrations 
    flux    = out.flux;                         % fluxes
    gasex   = out.gas;                          % gas exchange fluxes
    mix     = out.mix;                          % ocean mixing
    T       = out.T;                            % surface temperature
    RF      = out.rf;                           % radiative forcings
    t       = out.t;                            % time
    inp     = out.inp;                          % input parameters
    indx    = out.indx;                         % species/reservoir index
    tdep    = out.tdep;                         % time dependent forcings
    tau     = out.tau;                          % timescales
end