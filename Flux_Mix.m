%% ============================= EONS Model =============================== 
% Julia Horne, 2017

% Ocean mixing and sediment diffusion for all dissolved species. 

function mix = Flux_Mix(r,v)

spec = {'N2','RN','HNO3','H3PO4','DIC','TA','O2','CH4','FeO'};
bsd = {'s','d'}; bnz = {'n','z'}; 

for ii = 1:length(spec)
    for ib = 1:length(bsd)
       c.(bsd{ib}).(spec{ii}) = r.(bsd{ib}).(spec{ii})./v.oc.vol.(bsd{ib}); % mol/m3 
    end
% ocean overturning mixes dissolved species according to concentration
% gradient at an assumed rate (oc.mix.sd)
     mix.(spec{ii}).s = v.oc.mix.sd .* (c.d.(spec{ii}) - c.s.(spec{ii}));   % mol/yr; mixing in/out surface box
     mix.(spec{ii}).d = v.oc.mix.sd .* (c.s.(spec{ii}) - c.d.(spec{ii}));   % mol/yr; mixing in/out deep ocean
 
    for ib = 1:length(bsd) 
        b_oc = bsd{ib}; % ocean box
        b_se = bnz{ib}; % sediment box
% dissolved species diffuse through sediments into overlying ocean box, at 
% a rate controlled by the concentration gradient. Positive is into seds,
% negative is into ocean. 
    switch spec{ii}
        case 'FeO' % no diffusion for you!
        otherwise 
           c.(b_se).(spec{ii}) = r.(b_se).(spec{ii})./v.oc.vol.(b_se);      % mol/m3
           mix.diff.(spec{ii}).(b_se) = v.sed.diff.(spec{ii}) .* ((c.(b_oc).(spec{ii}) - c.(b_se).(spec{ii}))...
               ./ v.sed.difflength) .* v.sed.sa.(b_se);                     % mol/yr; diffusion of dissolved species into sediments
    end
    end
end
end
