%% ============================= EONS Model =============================== 
% Julia Horne, 2019

% Calculate CO2 and temperature dependent weathering modifiers.
% Weathering of silicate and carbonates is enhanced by higher temperatures
% and limited by substrate transport (ie. how accessible fresh rock is at 
% surface; Mills et al. 2011). We modify the weathering enhancements from
% COPSE (Bergman et al. 2004) that are due to pCO2, temperature, and plant
% colonization such that they are normalized to ~ 1 at the modern day.
% Seafloor weathering of basalts is here assumed to be proportional to the
% dissolved [CO2] in deep ocean/seds relative to modern day (Coogan +
% Gillis, 2013)

function [ws] = WeatheringSensitivities(t,T,CO2,dCO2,tdep,v)
refsz = size(t); 

sil = zeros(refsz); carb = zeros(refsz); sfw = zeros(refsz); dicsens = zeros(refsz); dCO2_PL = zeros(refsz);
CO2_PAL = zeros(refsz); Tsil = zeros(refsz); Tcarb = zeros(refsz); Tbasa = zeros(refsz); Wmod = zeros(refsz);

Econt = 20.5;                                               % kcal/mol; average activation energy for diopside (Brady + Carrol, 1994, Schott et al. 1981)
Ebasa = mean([11.7, 19]);                                   % kcal/mol; average activation energy for enstatite pyroxene + forsterite olivine (Brady + Carrol, 1994; Wogelius + Walther, 1992)
Ecarb = 32;                                                 % kJ/mol; activation energy for dolomite dissolution (Herman + White, 1985)
Rcal = (v.const.R./v.conv.Jcal).*1e-3;                      % kcal/Kmol; gas constant converted from joules
RkJ = v.const.R.*1e-3;                                      % kJ/Kmol; gas constant converted to kilojoules
tex = 0.5;                                                  % exponent (max 0.5, min 0.2) for CO2 sensitivity of terrestrial fluxes
oex = tex;                                                  % exponent for CO2 sensitivity of ocean flux

% temperature/pCO2 silicate/carbonate weathering modifiers adapted from COPSE (Bergman 2004) and Brady (1991) and Rushby et al (2018)
for ia = 1:length(t)
    CO2_PAL(ia) = CO2(ia)./v.atm.CO2pal;                    % PAL; fraction of modern CO2 levels (ie. a in COPSE)
    dCO2_PL(ia) = dCO2(ia)./v.oc.dCO2pl;                    % present level; fraction of modern z seds [CO2] 
    dicsens(ia) = dCO2_PL(ia).^oex;                         % sensitivity of seafloor weathering on [CO2]
    Tsil(ia) = exp((Econt/Rcal).*((1/289)-(1./T(ia))));     % Arrhenius equation for elevated dissolution of silicates (Brady + Carroll 1994)
    Tcarb(ia) = exp((Ecarb/RkJ).*((1/289)-(1./T(ia))));     % Arrhenius equation for elevated dissolution of carbonates 
    Tbasa(ia) = exp((Ebasa/Rcal).*((1/289)-(1./287.15)));   % Arrhenius equation for elevated dissolution of basalts assuming hydrothermal water is ~ 10C warmer than deep ocean
    % Endmember weathering modifier is either purely pCO2/temp controlled
    % (pre evolution of land plants, tdep.plant == 0) or is controlled by
    % plant fertilization and enhanced weathering at high pCO2,T
    % (tdep.plant == 1); treatment in COPSE/GEOCARB is a bit fudged to
    % prevent high pCO2 pre-plants, but in truth the modifier SHOULD
    % INCREASE over the plant evolution transition. So we just force the
    % pCO2 relationship to be stronger with plants by a factor of 2 (ests.
    % for plant chemical weathering enhancement is ~ 2-4 fold increase
    % across transition from minimal to heavily vegetated; Ibarra et al.
    % 2019- lower est.; Berner 2006- higher est.)
    
    Wmod(ia) = (1+tdep.plant(ia)) .* (CO2_PAL(ia).^tex); 	% weathering plant/CO2 sensitivity
    
    % temperature sensitivity is multiplied by weathering modifier. Wmod 
    % ref = 2 since we have plants today and CO2 = 1 PAL, normalizing this
    % variable.
    sil(ia) = Tsil(ia) .* (Wmod(ia) ./ 2); 
    carb(ia) = Tcarb(ia) .* (Wmod(ia) ./ 2);
    sfw(ia) = Tbasa(ia) .* dicsens(ia); 
end 

% weathering through time has a potential dynamic range of 1-10 (Mills et
% al. 2011, fig. 3)
b = 9; a = 1+b; 
ws.carb = (a .* carb) ./ (carb + b);                        % modern reference is < 2.4x maximum theoretical Phanerozoic rate
ws.sil  = (a .* sil) ./ (sil + b);   
ws.sfw  = sfw ./ (1 + sfw) ;                                % normalized to modern == 1

end