%% ============================= EONS Model =============================== 
% Julia Horne, 2022

% Sums all reservoirs of the same species across all ocean/sediment boxes 
% SO YOU DON'T HAVE TO. 
% inputs are the reservoir structure (rin == r) from the model output, and a
% list (cell!) of species names that you want to be totalled in all ocean+sediment
% reservoirs (ie. {'H3PO4'} == H3PO4 in s,d,n,z). 

function totres = TotalOceanReservoirs(rin,spes)
boxes = {'s','d','n','z'}; 

for ix = 1:length(spes)
    totres.(spes{ix}) = 0; % initialize total reservoir output
    switch spes{ix}
        case {'fixN','fixedN'} % total of RN + HNO3
            for ibx = 1:length(boxes)
                totres.(spes{ix}) = totres.(spes{ix}) + rin.(boxes{ibx}).RN + rin.(boxes{ibx}).HNO3; 
            end
        case 'FeO'
            bxz = {'s','d'}; 
            for ibz = 1:length(bxz)
                totres.(spes{ix}) = totres.(spes{ix}) + rin.(bxz{ibz}).(spes{ix}); 
            end
        otherwise % add up the individual reservoir fluxes
            for ibx = 1:length(boxes)
                totres.(spes{ix}) = totres.(spes{ix}) + rin.(boxes{ibx}).(spes{ix}); 
            end
    end
end

end
