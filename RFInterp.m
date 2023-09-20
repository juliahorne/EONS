%% Interpolate GHG radiative forcing based on mixing ratio
% Julia Horne, 2018

% Function calculates RF for a gas from the line-by-line transfer
% model provided by Brendan Byrne (2014). Reads in the gas abundance (mol
% gas / mol atm) from the EONS model, interpolates onto model output 
% to find a correlated RF value for those concentrations.

function [rad] = RFInterp(xGHG,gasname,bb)

% interpolate the model output into the SI data (bb.__)
% interp1(x,v,xq) where x = sample points, v = data vector, xq = query
% points (ie. what I want to interpolate into the vector)
logGHG = log10(xGHG);            % evaluate the log of the GH gas mixing ratio
xxGHG  = real(logGHG);           % ensure that there are no imaginary components!
rad = interp1(log10(bb.c.(gasname)),bb.rf.(gasname),xxGHG,'linear','extrap'); 

end
