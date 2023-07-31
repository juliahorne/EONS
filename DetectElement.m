%% determine if an element exists in a species, and count it's occurance 
% Julia Horne, 2021

% check for an element (EL) and it's amount in a species (SP) 
% returns a boolean t/f and (if present) the number of element in that
% species (ie. CO2 searching for element O would return [1,2]) 
% NOTE: this has limited utility beyond the species included in the EONS
% model, but should work on most simple chemical formulae. 

function [isELthere,mult] = DetectElement(SP,EL) 
    % now search for all instances of the element in species names 
    isELthere = contains(SP,EL); % returns boolean t/f
    if isELthere == 1 
         % first check that the element is not two characters long! make 
         % the element index reflect that by adding to the indexes
        if length(EL) == 2
            addon = 1; 
        else 
            addon = 0; 
        end
        % also check that there is not another element in the species that
        % contains that letter (ie. Ca and C) by checking if the trailing
        % value is a lowercase letter
        if length(strfind(SP,EL)) > 1                       % if there is more than one possible location for the element in the species name...
            xind = strfind(SP,EL); 
            for ix = 1:length(xind)
                if ischar(SP(xind(ix)+1))                   % check if the following value is a letter 
                    islo = upper(SP(xind(ix)+1))~=SP(xind(ix)+1); % check if it is lowercase
                    if islo == 1 
                        % if this is a lowercase letter, then we are looking at a different element... 
                    else % the next value is either a number or an uppercase letter, so this is the index we want
                        elind = xind(ix) + addon; 
                        mult = FindMult(SP,elind); 
                        return 
                    end
                end
            end
        end
        % check if there are multiple of the elements in this species
        % we typically have OH as batch multiplied in some species (ie.
        % B(OH)4) so make sure that O as well as H is multiplied by the
        % final number! Also check for parentheses, as a clue to where the
        % multiplier is!
        if contains(SP,'OH') && strcmp(EL,'O')              % if the species contains OH in this order, assume the next number after H is the multiplier (ie. (OH)_2)
            % only do this if the desired element is O!
            elind = strfind(SP,'H') + addon;
        elseif contains(SP,')') % check for any parentheses!
            par1ind = strfind(SP,'('); par2ind = strfind(SP,')'); % find the indexes of the two parentheses
            elind = strfind(SP,EL);                         % and the index of the element
            numafter = str2double(SP(elind+1));
            % if there is a number directly after the species, then consider that in the multiplier calculation
            if ~isnan(numafter) && elind > par1ind && elind < par2ind  
                spenum = str2double(SP(elind+1));           % the number in the parentheses after the species
                parnum = str2double(SP(par2ind+1));         % the number after the trailing parenthesis
                mult = spenum .* parnum; 
                return % exit to calling function!
            % if the species is in between those parentheses with no included number, then the index after the last parenthesis is the multiplier
            elseif elind > par1ind && elind < par2ind 
                elind = strfind(SP,')') + addon;
            % if the species lies outside of the parentheses, then ignore them!
            elseif elind < par1ind || elind > par2ind 
                elind = strfind(SP,EL) + addon; 
            end
        else % just find the index of the element
            elind = strfind(SP,EL) + addon; 
        end
        % the following character, if a number, is the multiplier
        mult = FindMult(SP,elind); 
    else % if not there, then return 0 multiplier
        mult = 0; 
    end

end

%% subfunction to find element multiplier

function mult = FindMult(SP,ind)
if strcmp(SP(ind),SP(end))                                  % make sure the index is not the end of the string
        mult = 1; % if there is no number, the multiplier is 1
else % proceed!
    mult = str2double(SP(ind+1));                           % make a multiplier out of the number
    if isnan(mult) % if it can't be made into a number, then it's a letter!
        mult = 1; 
    end
end
    

end
    