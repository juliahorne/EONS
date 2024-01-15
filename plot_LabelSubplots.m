%% Generate subplot labels, numerical or alphabetical
% Julia Horne, 2023

% based off of Sterling Baird's response:
% https://www.mathworks.com/matlabcentral/answers/435513-how-can-i-label-my-graphs-as-a-b-c-etc-in-subplot-matlab

% necessary input is axis handle (can be a structure)

% Optional inputs include location of labels (loc == 'tr' for top right,
% 'tl' for top left, 'br' for bottom right, and 'bl' for you get the idea)
% letters/numbers (alphnum = 'alpha','num','letter','number') and fontsize.
% Another optional input is the number of columns of subplots - this is
% really only important to designate if you are plotting a single column of
% subplots (col), otherwise the spacing works for 2-5 columns. 
% Defaults are top left alphabetical font size 9 without alignment.
% You can also opt for outside the figure (but this is only for top left,
% next to the yaxis label; loc = out). 

function plot_LabelSubplots(ax,loc,alphnum,fsz,col)
    if ~exist('loc','var')
        loc = 'tr'; 
    end
    if ~exist('alphnum','var')
        alphnum = 'alpha'; 
    end
    if ~exist('fsz','var')
        fsz = 9;
    end
    if ~exist('col','var')
        col = 2; 
    end
    % subtitle function is only available for versions > 2020b, so check
    % the version and use a different function if earlier
    v = version; 
    if contains(v,'R2020b')
        STITLE = 'Subtitle'; 
    else
        STITLE = 'Title';
    end
    % generate a list of labels first
    switch alphnum
        case {'alpha','letter'} 
            labs = ('a':'z')';
        case {'num','number'}
            labs = ('1':'9')''; 
        otherwise 
            error('Please designate either alphabetical (alpha, letter) or numerical (num, number) for label'); 
    end
    if col == 1
        colx = 0.5855;                         % for a single column of subplots, needs to be < 0.6
    else
        colx = 0.665; 
    end
    % loop through all axes in the input
    for ia = 1:length(ax)
        % get subtitle position = [x y height?] - always originate in the top center above the plot!
        ax(ia).(STITLE).Units = 'normalized';  % important! sometimes will switch to 'data'!!
        stpos = ax(ia).(STITLE).Position;   
        
        if contains(loc,'t')                   % determine top or bottom
            ypos = stpos(2) - 0.05; 
        else
            ypos = stpos(2) - 0.95; 
        end
        if contains(loc,'l')                   % determine left or right 
            xpos = stpos(1) - 0.45; 
        else
            xpos = stpos(1) + 0.45; 
        end
        if contains(loc,'o')                   % designates outside, top left above y label 
            ypos = stpos(2) - 0.05; 
            xpos = stpos(1) - colx;            % may need to be adjusted based on how one squeezes one's plots
        end
        ax(ia).(STITLE).Position(1) = xpos;     ax(ia).(STITLE).Position(2) = ypos;
        ax(ia).(STITLE).String = compose("(%s)",labs(ia));
        ax(ia).(STITLE).FontSize = fsz;         ax(ia).(STITLE).FontWeight = 'bold'; 
    end

end

