%% Calculate ocean concentrations in mol/m3!!
% Julia Horne, 2023

% inputs are reservoirs/concentrations and the x and v constant structures.

function c = VolumetricConcentrations(r,conc,v)
    bx = {'s','d','n','z'};                                         % evaluate all ocean/sed boxes
    for ie = 1:length(bx)
        Rsps = fieldnames(r.(bx{ie}));                              % all the reservoir species in this box
        Csps = fieldnames(conc.(bx{ie}));                           % all the concentrated species in this box
        sps  = unique([Rsps;Csps]);                                 % all unique species in this box
        for is = 1:length(sps)
            nm = sps{is};   
            if isfield(conc.(bx{ie}),nm)
                c.(bx{ie}).(nm) = conc.(bx{ie}).(nm) .* v.oc.rho;   % mol/m3; convert from mol/kg
            else
                c.(bx{ie}).(nm) = r.(bx{ie}).(nm) ./ v.oc.vol.(bx{ie}); % mol/m3
            end
        end
        % make a fixed N concentration
        c.(bx{ie}).fixN = c.(bx{ie}).RN + c.(bx{ie}).HNO3;  
    end

end
