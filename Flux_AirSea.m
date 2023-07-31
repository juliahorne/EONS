%% ============================= EONS Model =============================== 
% Julia Horne, 2018
%
% Air-sea gas exchange flux equations for gaseous species.
% NOTE: positive is into atmosphere and negative is into ocean 

function gasex = Flux_AirSea(r,conc,T,v)
T0                = 298.15;                                    % K; reference temperature for Henry's Law eqn 
spe               = {'N2','CO2','NH3','CH4','O2'};

for iz = 1:length(spe)
    p.(spe{iz})   = v.atm.P.*(r.a.(spe{iz})./v.atm.mol);       % Pa; partial pressure
    Kh.(spe{iz})  = v.atm.K0.(spe{iz}).*exp(v.atm.vH.(spe{iz})*(1./T-1/T0));% mol/m3Pa; Henry's gas-liq proportionality constant for temp T (kappa)
    switch spe{iz}
        case {'CO2','NH3'} % these are implicit species in the ocean
            c.s.(spe{iz}) = conc.s.(spe{iz}) .* v.oc.rho;      % mol/m3; concentration converted from mol/kg
        otherwise
            c.s.(spe{iz}) = r.s.(spe{iz}) ./ v.oc.vol.s ;      % mol/m3; concentration of X in stagnant boundary layer
    end
    c.a.(spe{iz}) = p.(spe{iz}) .* Kh.(spe{iz});               % mol/m3; concentration of X in atm
    D.(spe{iz})   = (v.atm.D.(spe{iz})*v.conv.spyr)/(100^2);   % m^2/yr; diffusion constant for gas X
    k.(spe{iz})   = D.(spe{iz})/v.oc.stag;                     % m/yr; gas-phase transfer velocity across stagnant boundary layer 'z'
    % air/sea gas exchange flux   
    gasex.(spe{iz}) = k.(spe{iz}) .* v.oc.sa .* (c.s.(spe{iz}) - c.a.(spe{iz})); % mol/yr
end

end
