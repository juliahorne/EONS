
%% ============================= EONS Model =============================== 
% Julia Horne, 2021

% Check that the model is conserving mass through time, particularly in H
% and O. The true change in tracked model elements SHOULD be closely
% related to the absolute tolerance (abstol in ops). We expect residual H
% and O to be misssing beyond tolerance because we don't explicitly track
% every single mole of H2O in the system. If the raw model output change in
% these two elements is in a 2:1 ratio, then we assume that untracked water
% is the cause. Any deviation from a perfect 2:1 ratio is cause for alarm. 

%% First, unpack the model run and sum up all the individual species across 
% all reservoirs
[t,r,conc,flux,gasex,mix,rf,T,~,tdep,~,~] = UnpackOutput(out);
% put reservoir and concentration output into the SumAllSpecies function to
% calculate the total moles of C,N,P,H,O,Fe,etc. in the model through time
[totalres, specres, list, roc, rsp] = SumAllSpecies(r,conc,inp.indx,v); 

% report total differences from all species
for itot = 1:length(list.elements)
    el = list.elements{itot}; % each of the major elements we are interested in
    % calculate the change
    totdif.(el) = totalres.(el)(end) - totalres.(el)(1);
    disp(['Change in tracked total ',el,' = ',sprintf('%0.3e',totdif.(el)),' mol']); 
end

% calculate change across all elements, plot up preliminary results
[change] = plot_MassConservation(t,totalres,totdif,list,v,'Raw'); 

%% Now compensate for O and H change from the dissociation of H2O 
ratio_OH = totdif.H./totdif.O ; 
disp ' '; % spacing!
disp(['Total H : O change across model run = ',num2str(ratio_OH)]); 
% generate "preliminary" total O and H reservoirs, changes
totalres.Hprelim = totalres.H; 
totalres.Oprelim = totalres.O; 
change.Hprelim   = change.H; 
change.Oprelim   = change.O; 
% then add this unbalanced mass to an implicit H2O reservoir, and add those
% H and O mols to the total reservoirs
roc.H2O = abs(change.O); % mol H2O
if (ratio_OH < 2.05 ) && (ratio_OH > 1.95) % reasonable range to assume that H2O is to blame
    if (totdif.O < 0) && (totdif.H < 0) % if O, H are lost, add the difference
        totalres.H = totalres.H + 2.*roc.H2O;
        totalres.O = totalres.O + roc.H2O; 
    elseif (totdif.O > 0) && (totdif.H > 0) % if O, H are gained, remove the difference
        totalres.H = totalres.H - 2.*roc.H2O; 
        totalres.O = totalres.O - roc.H2O;
    end
    ohlist = {'O','H'}; % recalculate true change in O and H
    disp ' '; % spacing for readability
    for itot = 1:length(ohlist)
        el = ohlist{itot}; 
        % recalculate the change between each timestep for O and H
        change.(el) = totalres.(el)-totalres.(el)(1); 
        % recalculate the total difference in the O, H reservoirs
        totdif.(el) = totalres.(el)(end) - totalres.(el)(1);
        % recalculate the ratio of change in the O, H reservoirs
        difratio.(el) = change.(el) ./ totalres.(el)(1) ; % the ratio of total change flux (X_t - X_t=1 / X_t=1)
        totdifratio.(el) = totdif.(el) ./ totalres.(el)(1); % ratio of total difference (X_end - X_t=1 / X_t=1)
        % and report!
        disp(['Change in Raw total ',el,'       = ',sprintf('%0.3e',totdif.(el)),' mol']); 
        disp(['Change ratio in Raw total ',el,' = ',sprintf('%0.3e',totdifratio.(el))]);
    end
    % replot the O, H, and all other major elements' evolution
    [change] = plot_MassConservation(t,totalres,totdif,list,v,'Corrected'); 
elseif isnan(ratio_OH) % dividing by zero!
    disp('O and H are unchanged'); 
else 
    % end the evaluation!
end


%% Subfunction : plotting mass changes

function [change] = plot_MassConservation(t,totalres,totdif,list,v,comment)

% plot change in total reservoirs over time (treated results)
figure(); clf; hold on; 
for itot = 1:length(list.elements)
    el = list.elements{itot}; 
    % calculate the change between each timestep
    change.(el) = totalres.(el)-totalres.(el)(1); 
    plot(t,change.(el),'-o','color',v.color.(el),'DisplayName',sprintf(['Delta ',el,' = ',sprintf('%0.3e',totdif.(el)),' mol'])); 
end
set(gca,'xscale','log'); 
legend('-DynamicLegend','location','southwest'); 
xlabel('Time (yr)'); ylabel('Change in total element reservoir (\Delta mol)'); 
% symlog('y'); % bilaterally symmetric log scale on y axis
title({'Conservation of Matter (X(t) - X(t=1))'; [comment,' output']});

% Evaluate the change over time in total reservoirs for preliminary results
% (change flux should evolve towards net 0)
figure(); clf; hold on; 
dt = t(2:end) - t(1:end-1); tavg = 0.5.*(t(2:end)+t(1:end-1));
for itot = 1:length(list.elements)
    el = list.elements{itot}; % each of the major elements we are interested in
    plot(tavg,(totalres.(el)(2:end)-totalres.(el)(1:end-1))./dt,'-o','color',v.color.(el),'DisplayName',sprintf(['d',el,'dt'])); 
end
xlabel('Time (yr)'); ylabel('Change flux in total element reservoir (mol/yr)'); 
set(gca,'yscale','log','xscale','log'); 
% symlog('y'); 
legend('-DynamicLegend','location','southwest'); 
title({'Change flux in total model reservoirs'; [comment,' output']}); 


% change in mass of species as portion of total initial mass
figure(); clf; hold on; 
for itot = 1:length(list.elements)
    el = list.elements{itot}; % each of the major elements we are interested in
    difratio.(el) = change.(el) ./ totalres.(el)(1) ; % the ratio of total change flux (X_t - X_t=1 / X_t=1)
    totdifratio.(el) = totdif.(el) ./ totalres.(el)(1); % ratio of total difference (X_end - X_t=1 / X_t=1)
    plot(t,difratio.(el),'-o','color',v.color.(el),'DisplayName',sprintf(['d',el,'/',el,'(t=1) = ',sprintf('%0.3e',totdifratio.(el))])); 
end

xlabel('Time (yr)'); ylabel('Ratio of change in total element reservoir'); 
set(gca,'yscale','lin','xscale','log'); 
legend('-DynamicLegend','location','northwest'); 
title({'Ratio of change in total model reservoirs (\DeltaX(t) / X(t=1) )'; [comment,' output']}); 

end



