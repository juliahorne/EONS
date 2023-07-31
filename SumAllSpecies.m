%% ============================= EONS Model =============================== 
% Julia Horne, 2021

% Mass conservation is SUPER IMPORTANT. Make sure you're doing it by
% counting all species and all elements in the model through time.

% This function is called by the 'MassConservationn' script, but can be run
% on its own of course; this function also calls the 'DetectElement'
% function, which again can be used on it's own. The main inputs here are
% the reservoir and ocean concentration structures (model output), the
% reservoir/species index list, and the constants and standards structure. 

% Optional input for debugging is the 'dbinp' input; literally anything 
% string entered here will activate the debug subfunction and return a
% detailed output as the function searches and sums all the species. 

function [totalres, specres, list, roc, rsp] = SumAllSpecies(r,conc,indx,v,dbinp)

% determine if optional debugging parameter has been designated
if ~exist('dbinp','var') % if dbinp doesn't exist assign a rando value
    dbinp = 666; 
end 

totalres = [];                                                            % empty structure for storing total reservoirs (for elements)
specres = [];                                                             % empty structure for storing total reservoirs (for species)
boxes = unique(indx.res);                                                 % all reservoirs in the model
ocboxes = {'s','d','n','z'};                                              % ocean boxes

% first calculate all the C, H, and N species reservoirs from their
% concentrations (if there are DIC/RN reservoirs, in case we are running a pared-down version of the model) 
alltrackedspecies = unique(indx.sp);                                      % all tracked species in the model
contentsDIC = contains(alltrackedspecies,'DIC'); 
contentsRN = contains(alltrackedspecies,'RN'); 
contentsP  = contains(alltrackedspecies,'H3PO4'); 


if any(contentsDIC) && any(contentsRN)                                    % checks if there are any instances of boolean T in an array
    if ~any(contentsP)                                                    % no P in the model
        elementlist = {'C','N','O','H'};                                  % major elements in the model, neglecting P and Fe
        totalres.P = 0;
    else
        elementlist = {'C','N','P','O','H','Fe'};                         % major elements in the model, including P and Fe
    end
    sps = {'CO2','HCO3','CO3','OH','H','NH4','NH3','BOH3','BOH4'};        % implicit ocean species
    for is = 1:length(sps)
        roc.(sps{is}) = 0;                                                % initialize total count
        for ib = 1:length(ocboxes)
            rsp.(ocboxes{ib}).(sps{is}) = v.oc.m.(ocboxes{ib}).*conc.(ocboxes{ib}).(sps{is}); % implicit reservoir for that box
            roc.(sps{is}) = roc.(sps{is}) + rsp.(ocboxes{ib}).(sps{is});% whole ocean implicit reservoirs
            debug_sumallspecies(dbinp,'res',sps{is},ocboxes{ib},0);     % returns message if debugging
        end
    end
    allmodelspecies = [alltrackedspecies, reshape(fieldnames(roc),1,9)] ;
elseif any(contentsDIC) && ~any(contentsRN) && ~any(contentsP)            % no RN, Fe, or P reservoirs, only DIC
    sps = {'CO2','HCO3','CO3','OH','H','BOH3','BOH4'};                    % implicit ocean species
    elementlist = {'C','O','H'};                                          % major elements in the model
    totalres.P = 0; totalres.N = 0; 
    for is = 1:length(sps)
        roc.(sps{is}) = 0;                                                % initialize total count
        for ib = 1:length(ocboxes)
            rsp.(ocboxes{ib}).(sps{is}) = v.oc.m.(ocboxes{ib}).*conc.(ocboxes{ib}).(sps{is}); % implicit reservoir for that box
            roc.(sps{is}) = roc.(sps{is}) + rsp.(ocboxes{ib}).(sps{is});% whole ocean implicit reservoirs
            debug_sumallspecies(dbinp,'res',sps{is},ocboxes{ib},0);     % returns message if debugging
        end
    end
    allmodelspecies = [alltrackedspecies, reshape(fieldnames(roc),1,7)] ;
else % no DIC, RN reservoirs, so just use the tracked species in the model (ie. no speciation is happening, no implicit ocean reservoirs)
    allmodelspecies = alltrackedspecies;                                  % all species in the model
    for is = 1:length(allmodelspecies)
        xname = allmodelspecies{is};                                      % the species
        elist{is} = xname(1);                                             % add the first letter of that species name to the list of elements
    end
    elementlist = unique(elist);                                          % only use the unique elements in the list!
end
specieslist = unique(allmodelspecies);                                    % all unique species in the model

%% there are some exceptions to be made for reservoirs that are tracked but more specifically speciated in the implicit ocean reservoirs
% remove these to make the DetectElement test experience fewer redundancies!

% make a special case for RN == NH4 + NH3 and DIC == HCO3 + CO3 + CO2
specieslist(strcmp(specieslist,'RN')) = [];                               % remove RN from the list of all species for the next evaluation
specieslist(strcmp(specieslist,'DIC')) = [];                              % remove DIC from the list of all species for the next evaluation
% create one ocean biomass reservoir
specieslist(strcmp(specieslist,'LB')) = [];                               % remove LB from the list of all species for the next evaluation


%% evaluate all species in the model across all boxes, summing them into one species reservoir

for isp = 1:length(specieslist)
    specname = specieslist{isp} ;
    ocsplist = fieldnames(roc);                                           % all species in the implicit ocean reservoirs
    specres.(specname) = 0;
    % first, check all of the explicit reservoirs for the species (all species reservoirs that have actual ODEs)
    for ib = 1:length(boxes)
       spslist = fieldnames(r.(boxes{ib}));                               % all species in that tracked reservoir
       % check that this reservoir contains the species (boolean t/f)
       isXinExplicitList = strcmp(spslist,specname);                      % is the species in the explicit reservoirs?
       for ireslst = 1:length(isXinExplicitList)
           if isXinExplicitList(ireslst) == 1                             % take that reservoir and add it to the total explicit reservoir
              specres.(specname) = specres.(specname) + r.(boxes{ib}).(specname);
              debug_sumallspecies(dbinp,'res',specname,boxes{ib},0);      % returns message if debugging
           else % move on
           end
       end
    end
    % next, check all of the implicit ocean reservoirs for the species
    % (species that aren't actually tracked, but derived from the speciation of explicit species in Flux_Spec)
    isXinImplicitList = strcmp(ocsplist,specname);                        % is the species in the implicit reservoirs?
    for ioclst = 1:length(isXinImplicitList)
        if isXinImplicitList(ioclst) == 1                                 % take that total ocean implicit reservoir and add it to the total 
            specres.(specname) = specres.(specname) + roc.(specname); 
            % no debugging for this line; already noted in the roc structure generation step
        else % move on
        end   
    end
end

%% to avoid any confusion and double-adding, explicitly create an organic carbon reservoir
if any(contains(specieslist,'OC'))                                        % if there is any OC in the model...
    specres.POC = 0;
    for iboc = 1:length(ocboxes)                                        % make a POC reservoir total that excludes the continental/u sed (ocean + reactive seds only)
       specres.POC = specres.POC + r.(ocboxes{iboc}).OC; 
    end
    if isfield(r.s,'LB')                                                  % add the LB reservoir to the total model organic carbon (buried and particulate)
        specres.OC = specres.OC + r.s.LB ;                                % total org carbon == ocean + sed + continental + live Org C
        specres.POC = specres.POC + r.s.LB;                               % ocean and live Org C
    end
    % things get a bit confusing for organic matter since we split into
    % organic P and N AFTER DEATH only, so add in the LB reservoir 
    nps = {'N','P'}; 
    for in = 1:length(nps)
        if any(contains(specieslist,['O',nps{in}]))                       % check if there is already a cont/sed ON or OP reservoir
            specres.(['O',nps{in}]) = specres.(['O',nps{in}]) + (r.s.LB./v.const.(['C',nps{in},'ratio']));
        else % if there's no sed or continent ON/OP, then the organic fraction of these species are from POC only
            specres.(['O',nps{in}]) = specres.POC ./v.const.(['C',nps{in},'ratio']);
            specieslist{end+1} = ['O',nps{in}];                           % add it to the species list!
        end
    end
end
   
%% Make sure that true stoichiometry is honored in the count!
stoilistall = specieslist;                                                % duplicate of the list of all model species
for ixs = 1:length(specieslist) 
    spe = specieslist{ixs};                                               % iterate through each species in the model 
    [specres,stoi] = continental_stoi(specres,spe,dbinp); 
    % with organic matter reservoirs all assigned to true stoichiometries,
    % we make a new specieslist for the next loop and ignore the
    % redundancies (OC, OP, ON)
    switch spe
        case {'OP','ON','OC','CP','SP'}
            stoilistall(strcmp(stoilistall,spe)) = {stoi};                % add true stoichiometry of organic components to the list
        otherwise % move on
    end
end
% the replacement of ON and OP with NH3 and H3PO4 leads to redundancies in this list. REMOVE!
stoilist = unique(stoilistall);                                           

%% then evaluate all species for each element, and add those to the total element reservoir
for iem = 1:length(elementlist)
    elem = elementlist{iem};                                              % which element in the model are we summing?
    totalres.(elem) = 0;
   for isx = 1:length(stoilist)
       spes = stoilist{isx};                                              % one species in the reservoir at a time
       ospes = specres.(spes);                                            % the total reservoir of the species 
       switch spes
           case 'FeIIP' % search for this species' ACTUAL formula, which cannot be the structure name (because parentheses)
               spsnm = 'Fe3(PO4)2';  
           otherwise
               spsnm = spes; 
       end
        % now search for all instances of the element in species names 
        [isXthere,mult] = DetectElement(spsnm,elem);                      % returns boolean t/f if the element is in the species, and how much of it is in there
        if isXthere == 1  
            if strcmp(elem,'H') && strcmp(spes,'HZ')                      % ignore the bullshit haze reservoir
                mult = 0; 
            elseif strcmp(spes,'FeIIP')                                   % measured in mols P, so divide by 2!
                mult = mult./2; 
            end
            % add that reservoir to the total
            totalres.(elem) = totalres.(elem) + mult.*ospes; 
            debug_sumallspecies(dbinp,'elems',elem,spsnm,mult);           % returns message if debugging
        else % move on!
            debug_sumallspecies(dbinp,'elems',elem,spsnm,mult);
        end
   end        
end


%% pack up lists for output
list.species = unique(stoilist); 
list.elements = elementlist; 

end

%% Subfunction : assign true stoichiometry to continental reservoirs for organic matter 
function [spres,stoi] = continental_stoi(specres,spe,dbinp)
spres = specres;                                                          % copy the specres structure to a new structure for posterity
% first assign a true stoiciometry to the inputted reservoir name 
switch spe
    case {'ON','OP'} % in this case, NH3 and H3PO4 already exist as inorganic reservoirs, so we don't start with a 0 stoi reservoir                                                      
        switch spe
            case 'ON' 
                stoi = 'NH3';
            case 'OP'
                stoi = 'H3PO4';
        end
        reslst = 'oc + u + c'; 
    case 'OC' % start with 0 stoi res because we don't track CH2O as a model species
        stoi = 'CH2O';
        spres.(stoi) = 0;
        reslst = 'oc + u + c';
    case {'SP','CP'} % inorganic P reservoirs (silicate + carbonate) == PO4, FePO4
            stoi = 'PO4';
        if isfield(spres,stoi) % if this is already been initialized, ignore!
        else
            spres.(stoi) = 0;
        end
        reslst = 'u + c';
    otherwise 
        stoi = spe;  
        return                                                            % return to parent function
end
% now make a new specres entry (or add to existing ones..) for those 
% reservoirs under the true stoichiometry 
spres.(stoi) = spres.(stoi) + spres.(spe);
debug_sumallspecies(dbinp,'res',stoi,reslst,0);                     % returns message if debugging
end



%% Subfunction : debugging output messages
function debug_sumallspecies(dbg,desig,e,s,m)
if strcmpi(dbg,'debug')                                                   % execute if flag is true
    if strcmpi(desig,'elems')                                             % detail how the individual elements are being summed in each species
        message1 = [e,' in ',s,' ',num2str(m),' times'];                  % display element, species, and multiple times element is in species
        switch s 
        case 'TA'
            switch e
            case {'C','N','O','H'}
                message2 = '** already counted in implicit species'; 
                disp([message1,' ', message2]);
            otherwise
                disp(message1); 
            end
        case 'OC'
            switch e
                case {'P','N'}
                    if strcmp(e,'N')
                        message2 = '** already counted in ON' ;
                    elseif strcmp(e,'P')
                        message2 = '** already counted in OP' ;
                    end
                    disp([message1,' ', message2]); 
                otherwise
                    disp(message1); 
            end
        case {'ON','OP'}
            switch e
                case {'H','O'} 
                    message2 = '** already counted in OC'; 
                    disp([message1,' ', message2])
                otherwise 
                    disp(message1); 
            end
            otherwise % print normal message
                disp(message1); 
        end
        
    elseif strcmpi(desig,'res')                                           % detail how we count each species in each reservoir of the model
        disp([e,' in ',s,' reservoir']);                                  % display reservoir and the species within it
    else 
    end
else 
end

end
