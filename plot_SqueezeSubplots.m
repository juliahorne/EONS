%% Squeeze subplots together so close it hurts
% Julia Horne, 2022

% Inputs are axes handles (axs = ax(1:n) == subplot(n,x,y)), and 
% order should designate whether or not the squeeze is vertical or
% horizontal (ord = 'v','h'). Also need to know if plotting only vertical
% or only horizontal! If you have a, say, 2x2 figure, the position of the
% second subplot will be different! Input 'pos' should be the horizontal
% subplot length if squeezing vertically. If squeezing horizontally, count
% will be the same (ie. for 2x2, 'pos' = 2 for 'v' and = 1 for 'h').

% Optional inputs (but a good idea) are spacing (spacing, duh!), yes/no
% width/height change (adwh) and presence of axis labels x or y; default 
% spacing = 'comfort', so we we leave more room in between subplots, which
% is the same for 'y' = axlbl, for some extra spacing for a label. Default
% axlbl if not entered is 'n' (so tighter spacing with no axis labels).
% Default for spacing is 'comfort', medium-tight spacing (options include
% 'comfort','tight','none', or 'normal'). 

% You can also request 'widen', which is NOT widening intraplot spacing but
% widening the plot itself, via x position. BEWARE: this doesn't take y 
% axis labels into account, so be prepared to only have space for a single
% line of text outside of the axis ticks if fontsize is > 9. This is really
% only useful if you're otherwise squeezing vertically.

function plot_SqueezeSubplots(axs,ord,pos,spacing,axlbl,widen)
if ~exist('spacing','var')
    spacing = 'normal';
end
if ~exist('axlbl','var')
    axlbl = 'n'; 
end
if ~exist('widen','var')
    widen = 'n'; 
end

% determine which axes dimensions to manipulate (Position = [x, y, width, height])
switch ord
    case 'v' % vertical dimensions
        xy = 2; wdht = 4; % y and height
        roomy = 1/9.5; tight = 1/6; % roomy and tight fractions
        axis = 'x'; axnum = length(axs); % the axis labels to remove and the one to remain
    case 'h' % horizontal dimensions
        xy = 1; wdht = 3; % x and width
        roomy = -1/20; tight = -1/5; % roomy and tight fractions
        axis = 'y'; axnum = 1; % the axis labels to remove and the one to remain
    otherwise
        error('Second input ORDER must be either v (vertical), or h (horizontal) squeeze direction');
end
switch spacing
    case 'comfort' % more space!
        frac = roomy; 
    case 'tight' % very little space!
        frac = tight;
    case 'none' % obviously this means no space between subplots!
        frac = 1/3.5; 
    case 'normal' % no change!
        frac = 0; 
    otherwise 
        error('Fourth input SPACING must be comfort (hella space), tight (self-explanatory), normal (no change), or none (no gaps between plots)');
end
switch ord
    case 'v'
        hfrac = 1; vfrac = abs(frac);
    otherwise 
        hfrac = abs(frac); vfrac = 1;
end

    % if it is chosen to widen the subplots, then make that adjustments 
    % before the squeezing step (moving left, up, and widening)
    if strcmp(widen,'y')
       initialleftmostpos = axs(1).Position(1);
       for ii = 1:length(axs)
           % if there are 2 columns, change the left adjustment to a right adjustment
           if axs(ii).Position(1) > initialleftmostpos % triggers when the x position is larger!
              axs(ii).Position = axs(ii).Position + [0.025 0.0275 0.1225/pos 0]; 
           else
              axs(ii).Position = axs(ii).Position + [-0.035 0.0275 0.1225/pos 0]; 
           end
       end
    end

    % displacement between subplots is uniform, so adjust all in a
    % consistent manner in x/y space. Numbering changes if plots are
    % vertically or horizontally aligned! 
    dist = axs(1).Position(xy) - axs(1+pos).Position(xy); % intraplot distance
    
    for ii = 1:length(axs)
        axs(ii).Position(wdht) = axs(ii).Position(wdht) + frac*dist; 
        if strcmp(axlbl,'n') % remove axis ticks if desired
            if ii ~= axnum
                axs(ii).([upper(axis),'TickLabels']) = {[]};
            end
        end
        % then adjust all plots slightly down/to the right by a portion of
        % the increased height/width to center
        axs(ii).Position(xy) = axs(ii).Position(xy) - (1/2).*frac*dist; 
        % AND in the case that the number of subplots is > 6, slightly
        % decrease font size for easier reading
        if length(axs) > 6
            axs(ii).FontSize = axs(ii).FontSize - 2; 
        end
    end
    

end