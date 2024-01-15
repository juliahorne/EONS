%% ============================= EONS Model =============================== 
% Julia Horne, 2022

% Sums all fluxes of the same type across all ocean/sediment reservoirs SO
% YOU DON'T HAVE TO. 
% Inputs are the flux structure (infx == flux) from the model output, and a
% list (cell!) of flux names that you want to be totalled in all ocean+sediment
% reservoirs (ie. {'denit'} == denitrification of HNO3 in s,d,n,z). 
% EXAMPLE USE:
% totflx = TotalOceanFluxes(flux,{'denit','metha','mtrophy','ammon'});

function totflx = TotalOceanFluxes(infx,fxs)
boxes = {'s','d','n','z'}; 
for ix = 1:length(fxs)
    switch fxs{ix}
        case {'burial','sed','reduce'}
            totflx.(fxs{ix}) = []; % initialize as a structure
                spes = fieldnames(infx.(fxs{ix})); bnz = {'n','z'} ;
                for isp = 1:length(spes)
                    totflx.(fxs{ix}).(spes{isp}) = 0;
                   for ibz = 1:length(bnz)
                       totflx.(fxs{ix}).(spes{isp}) = totflx.(fxs{ix}).(spes{isp}) + infx.(fxs{ix}).(spes{isp}).(bnz{ibz}); 
                   end
                end
        otherwise
        totflx.(fxs{ix}) = 0;       % initialize as a value
        switch fxs{ix}
            case {'revweather','scav','sorb',}
                bnz = {'n','z'};
                for ibz = 1:length(bnz)
                    totflx.(fxs{ix}) = totflx.(fxs{ix}) + infx.(fxs{ix}).(bnz{ibz}); 
                end
            case {'forg','bureff'} % organic C burial fraction or burial OC efficiency
                spe = {'OC','CaCO3'}; bur = [];
                for ie = 1:length(spe)
                    bur.(spe{ie}) = infx.burial.(spe{ie}).n + infx.burial.(spe{ie}).z;
                end
                if strcmp(fxs{ix},'forg')
                    totflx.(fxs{ix}) = bur.OC ./ (bur.OC + bur.CaCO3); 
                else
                    totflx.(fxs{ix}) = bur.OC ./infx.prod.CO2;
                end
            otherwise
            for ibx = 1:length(boxes)
                switch fxs{ix} 
                    case 'axrm'
                        totflx.(fxs{ix}) = totflx.(fxs{ix}) + infx.denit.OC.(boxes{ibx}) + infx.metha.OC.(boxes{ibx}); 
                    case 'oxrm'
                        totflx.(fxs{ix}) = totflx.(fxs{ix}) + infx.ammon.OC.(boxes{ibx}); 
                    case 'denitN'
                        totflx.denit = totflx.(fxs{ix}) + infx.denit.HNO3.(boxes{ibx}); 
                    case {'ammon','metha','denitC'} 
                        if strcmp(fxs{ix},'denitC')
                            fxsX{ix} = fxs{ix}(1:end-1);
                        else
                            fxsX{ix} = fxs{ix}; 
                        end
                        totflx.(fxs{ix}) = totflx.(fxs{ix}) + infx.(fxsX{ix}).OC.(boxes{ibx}); 
                    otherwise % add up the individual reservoir fluxes
                        totflx.(fxs{ix}) = totflx.(fxs{ix}) + infx.(fxs{ix}).(boxes{ibx}); 
                end
            end
        end
    end
end

end
