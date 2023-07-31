%% Set a new automatic color order for a figure
% Julia Horne, 2022

% The number of auto colors produced by MATLAB is fixed, and that's very
% annoying, since plotting more than 7 lines will result in repeated colors
% if not otherwise specified. This function allows me to still autoplot
% without having to assign individual colors to each reservoir or flux,
% because that's tedious as shiiit.
% Inputs are the current figure/axis handle, the number of iterations to
% repeat (based on the 7 original random Matlab colors, producing 7*xcolors
% color triplets), and an option order that signals if the user wants a
% gradient, random, or a rainbow colormap! Personally, I recommend rainbow.

function newOrd = plot_ColorOrder(xgca,xcolors,order)
if ~exist('order','var') % an optional input
    order = 'random';                                       % generate a random colororder if no order is given
end
cOrd = xgca.ColorOrder;                                     % original color order from the current figure handle (blue,red,yellow,purple,green,cyan,magenta)
cbow = [cOrd(2,:);cOrd(7,:);cOrd(4,:);cOrd(1,:);cOrd(6,:);cOrd(5,:);cOrd(3,:)];  % reorder into a red->yellow rainbow
switch order                                                % what kind of colormap is generated is based on user input, if given
    case 'random'                                           % the modifications to each color triplet are randomized
        newOrd = cOrd;                                      % initialize the new color order matrix
        for inx = 1:xcolors
            newcols{inx} = rand(size(cOrd(1,:))); 
        end
        % randomly assign some of these modifiers to be negatives 
        % leads to many generated 0s, so we don't use this until generating a
        % significantly large color order size (wherein the colors get quite light
        % as it is)
        if xcolors > 6
            for inc = 4:length(newcols) 
                X = randperm(3);                            % random indices for the triplet values
                newcols{inc}(X(1:3)) = -(newcols{inc}(X(1:3))./2);% assign a small negative to some elements at the randomly generated 3 slot
            end
        end
        % now amend the color order matirx with my new modifiers
        for id = 1:length(newcols)
            for ico = 1:7
                newOrd(end+1,:) = cOrd(ico,:) + newcols{id};% should be 0 < x < 1
                for in = 1:length(newOrd(end,:))
                    xn = newOrd(end,in);                    % check each number isn't less than 0 or more than 1
                    if xn > 1
                        xn = xn./2;                         % divide by 2
                    elseif xn < 0
                        xn = xn - xn;                       % just make it 0
                    end
                    newOrd(end,in) = xn;                    % and replace
                end
            end
        end
        % the colors are randomly altered, but mix them up to ensure that any sampling is as varied as possible
        nm = length(newOrd(:,1)); Y = randperm(nm);        % the indexes of each generated triplet
        oldOrd = newOrd; 
        for iy = 1:nm  
            newOrd(iy,:) = oldOrd(Y(iy),:);                % move one triplet to another's position
        end
% this is amended from the technique shown in 'mycolortable.m' by DGM (https://www.mathworks.com/matlabcentral/answers/1622300-how-to-create-rainbow-colormap-with-violet)
    case 'rainbow'                                         % reorder, then expand between the matlab colors to make a rainbow
        na = size(cbow,1);
        newOrd = interp1(linspace(0,1,na),cbow,linspace(0,1,xcolors*length(cbow)));% interpolate xcolor times!
    case 'gradient'                                         % again, make a rainbow but make the colors fade from white
        cnew = zeros(length(cbow).*2,3); 
        for inc = 1:length(cnew)
            if rem(inc,2) == 1                              % a remainder of 0 indicates an even number
                cnew(inc,:) = [1 1 1];                      % add white between each color
            else
                cnew(inc,:) = cbow(inc./2,:);               % add the cbow color
            end
        end
        na = size(cnew,1);
        newOrd = interp1(linspace(0,1,na),cnew,linspace(0,1,xcolors*length(cbow)));% interpolate xcolor times!
        % repeats because of the interspersed white lines, so I need to
        % improve!
end   


%% plot up how it works!
% x = 1:10;
% for inx = 1:length(newOrd)
%    plot(x,x+inx,'Color',newOrd(inx,:)); hold on;  
% end

end