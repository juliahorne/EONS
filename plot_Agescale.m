%% Change from a timescale to an agescale (in gigannum or megannum)
% Julia Horne, 2022

% inputs are an axis handle (ax) and age units (either 'ga' or 'ma'); don't 
% make me do any actual work!

function plot_Agescale(ax,unit)
   xlabel(ax,['Age (',upper(unit(1)),'a)']);    % label it!
   switch unit
       case 'ga'    
           age = 4 - ax.XTick./1e9;           % calculate age from time in model run
       case 'ma'
           age = 4e3 - ax.XTick./1e6; 
       otherwise
           error('Unit must be either ga = gigannum or ma = megannum!'); 
   end
   % translate the ages into labels!
   for ix = 1:length(age)
      ax.XTickLabels{ix} = num2str(age(ix));  
   end

end