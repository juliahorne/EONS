%% ============================= EONS Model =============================== 
% Julia Horne, 2019

% Takes in the atmosphere reservoir structure, finds the partial pressure
% of GHGs and calculates an equilibrium temperature from the resulting 
% radiative forcings and the evolving solar constant. 

function [rf,Fs,Teq] = Flux_Temp(t,ra,inp,v)
gg     = {'NH3','CO2','CH4'};
nt     = t;

%% Calculate solar flux
% sun has increased in luminosity by 25-30% over 4.5e9 years; use 
% Caldeira + Kasting (1992) equation for solar luminosity increase. In
% their function, t = years FROM present (ie. negative time, hence
% subtracting the full duration of model run!)
if isfield(inp,'treset')
    nt = t + inp.treset;                     % use real model time if reset has occurred
end
rf.sc  = (1 - (0.38.*(nt - 4e9) ./4.55e9)).^(-1) .* v.S.Pref; % W/m2; solar luminosity 

Fs = (rf.sc./4).*(1-v.ea.alb);               % W/m2; solar energy influx

%% Calculate radiative forcings from greenhouse gases (gg)
% NOTE: Byrne + Goldblatt (2014) model output file needs to be on your path!
rf.t = 0;
for ig = 1:length(gg) 
    mr.(gg{ig}) = ra.(gg{ig}) ./ v.atm.mol;  % mixing ratio
    [rf.(gg{ig})] = RFInterp(mr.(gg{ig}),gg{ig},v.bb); 
    if rf.(gg{ig}) < 0 % at extremely low pNH3, the radiative forcing may go negative!
        rf.(gg{ig}) = zeros(size(rf.(gg{ig})));
    end
    rf.t = rf.t + rf.(gg{ig}); % add the radiative forcings together to get a total rf
end


%% Calculate the equilibrium temperature from all combined GHGs
Teq.t = (((((rf.sc/4).*(1-v.ea.alb))+(v.const.aRF.*rf.t.^v.const.bRF))./v.const.bol).^(1/4));

% temperature equilibria depending on each individual GHG
for ite = 1:length(gg)
    Teq.(gg{ite}) = (((((rf.sc/4).*(1-v.ea.alb))+(v.const.aRF.*rf.(gg{ite}).^v.const.bRF))./v.const.bol).^(1/4));
end 

end

