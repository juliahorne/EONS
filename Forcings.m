%% ============================= EONS Model =============================== 
% Julia Horne, 2022

% All time-sensitive forcings in the model are calculated here, except for
% matle influx (that's in the Flux_BGC function).

% Photosythesis evoles sometime in the Archean, here assumed to be at the
% start of the Eoarchean (3.6-3.5 Ga, ~500 million years ito the model run)

% Body size and oxygen are positively correlated (Alexander 1971, Payne et
% al 2011). We assume that the Archean starts with smaller organisms that
% sink out of the ocean at a slower rate. We further assume that at a
% certain oxygen level, arbitrarily chosen to be aroudn 1e-3 PAL, larger
% organisms can evolve. Thus at this oxygen level, OM sinking increases
% to the modern rate - this is implemeted as a "sink factor" that we
% multiply all of our ocean (s,d) sinking fluxes by. Modern rate is 1, 
% while a low oxygen rate is 5 (ie. modern sinking is 5x faster). We assume
% that large organisms, regardless of oxygen level, only appear after 4 
% billion years (Cambrian explosion!)

% There is also some evidence of enhanced cotinetal weathering at the
% Archean-Proterozoic boundary, likely via microbially-facilitated 
% dissolution processes, particularly with respect to silicate and 
% carbonate phosphate rocks (Bayon et al. 2022). We simulate this by a
% weathering timescale enhancement for inorganic P materials that
% evolves to the modern rate of P-weathering from these reservoirs
% after the colonization of land by lichen/fungi (>4 billion years)

% A fixed terrestrial net productivity flux turns on stepwise with fungi
% evolving (1/10th of modern land productivity) and then with the evolution
% of vascular plants. 

function tdep = Forcings(t,v)

tdep.sink = zeros(size(t)); tdep.fungi = zeros(size(t)); tdep.plant = zeros(size(t));
tdep.photo = zeros(size(t)); tdep.terr = zeros(size(t)); 

% Some factors evolve from 0-1 (ie. photosythesis is not operating, then
% is allowed to evolve with nutriets) versus some that are changing from
% factors of 2,10 to 1 (ie. timescales for weathering take longer without
% selective P mining)
for ip = 1:length(t)
    if t(ip) <= v.td.initphoto                     % all inactive
        tdep.photo(ip)    = 0;                      % photosythesis
        tdep.fungi(ip)    = 1/10;                   % land colonization by fungi 
        tdep.sink(ip)     = 1/5;                    % timed body size growth
        tdep.plant(ip)    = 0;                      % vascular land plant colonization
        tdep.terr(ip)     = 0;                      % terrestrial productivity
    elseif t(ip) > v.td.initphoto && t(ip) <= v.td.initfungi    % activate photosythesis
        tdep.photo(ip)    = 1 ./ (1 + 1e12.*exp(-v.td.photok .* (t(ip) - v.td.initphoto)));
        tdep.fungi(ip)    = 1/10;
        tdep.sink(ip)     = 1/5; 
        tdep.plant(ip)    = 0; 
        tdep.terr(ip)     = 0; 
    elseif t(ip) > v.td.initfungi && t(ip) <= v.td.initbody     % activate fungal land expansion and partial terrestrial biosphere
        tdep.photo(ip)    = 1 ./ (1 + 1e12.*exp(-v.td.photok .* (t(ip) - v.td.initphoto))); 
        tdep.fungi(ip)    = 1 ./ (1 + 9.*exp(-v.td.fungik .* (t(ip) - v.td.initfungi)));
        tdep.sink(ip)     = 1/5;
        tdep.plant(ip)    = 0;
        tdep.terr(ip)     = (1 ./ (1 + 1e12.*exp(-v.td.terrk .* (t(ip) - v.td.initfungi)))).*0.1; 
    elseif t(ip) > v.td.initbody && t(ip) <= v.td.initplant     % activate body size growth
        tdep.photo(ip)    = 1 ./ (1 + 1e12.*exp(-v.td.photok .* (t(ip) - v.td.initphoto)));
        tdep.fungi(ip)    = 1 ./ (1 + 9.*exp(-v.td.fungik .* (t(ip) - v.td.initfungi)));
        tdep.sink(ip)     = 1 ./ (1 + 4.*exp(-v.td.bodyk .* (t(ip) - v.td.initbody)));
        tdep.plant(ip)    = 0;
        tdep.terr(ip)     = (1 ./ (1 + 1e12.*exp(-v.td.terrk .* (t(ip) - v.td.initfungi)))).*0.1; 
    elseif t(ip) > v.td.initplant                         % activate vascular land plant evolution and complete terrestrial biosphere
        tdep.photo(ip)    = 1 ./ (1 + 1e12.*exp(-v.td.photok .* (t(ip) - v.td.initphoto)));
        tdep.fungi(ip)    = 1 ./ (1 + 9.*exp(-v.td.fungik .* (t(ip) - v.td.initfungi)));
        tdep.sink(ip)     = 1 ./ (1 + 4.*exp(-v.td.bodyk .* (t(ip) - v.td.initbody)));
        tdep.plant(ip)    = 1 ./ (1 + 1e12.*exp(-v.td.plantk .* (t(ip) - v.td.initplant)));
        tdep.terr(ip)     = (1 ./ (1 + 1e12.*exp(-v.td.terrk .* (t(ip) - v.td.initfungi)))).*0.1...
            + (tdep.plant(ip) .* 0.9); 
    end

end

end

