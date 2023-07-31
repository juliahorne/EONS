%% Make a bunch of subpots with boxes ON
% Julia Horne, 2022

% inputs are the number (num) of subplots desired, and the direction of
% their arragnement (dim = 'h' or 'v', horizontal or vertical OR 2 values 
% designating y,x lengths). Optional inputs are subplot height (ht) and
% merged subplot numbers (ie. nmrg = [1 2]). If nmerge is designated, then
% those subplots will be merged into one big subplot, and you need to add
% one extra plot to your count (num = num+1). Can enter more than one set
% of subplot positions to merge. 
% 
% EXAMPLE USE:
%   [f,ax] = GenerateSubplots(5,'v');       == 5 vertical subplots
%   [f,ax] = GenerateSubplots(5,'h');       == 5 horizontal subplots
%   [f,ax] = GenerateSubplots(4,[2 2]);     == 4 subplots, 2 by 2
%   [f,ax] = GenerateSubplots(5,'v',[1 2]); == 4 vertical subplots, first subplot 2x tall
%   [f,ax] = GenerateSubplots(6,[3 2],{[1 2],[3 4]}); == 3 x 2 subplots, first and second subplot 2x wide (ie. 4 subplots total)

function [f,ax] = plot_GenerateSubplots(num,dim,nmrg,ht)
    if ~exist('nmrg','var')
        nmrg = 'x';
    end
    if ~exist('ht','var')
       ht = 'x';  
    end
    if ischar(dim)  % either vertical or horizontal, single column
        switch dim
            case 'v'
                dims = [num 1]; 
            otherwise
                dims = [1 num];
        end
    else 
        dims = dim; % inputs should be in the order [# vertical, # horizontal]
    end
    f = figure(); hold on;
    if ischar(nmrg) % continue as normal
        for i = 1:num
            ax(i) = MakeSubplotDims(dims,i);
            ChangeHeight(ax(i),ht); 
        end
    else
        for i = 1:num-1 % count one less axis!
            if iscell(nmrg)
                for im = 1:length(nmrg)
                    mrgn = nmrg{im}; % do the first group merge in the list
                    ax(i) = MergeSubplots(mrgn,dims,ht,i);
                end
            else 
                ax(i) = MergeSubplots(nmrg,dims,ht,i); 
            end
        end 
    end
end

%% subfunction : assign subplot dimensions
function axs = MakeSubplotDims(dim,i)
    axs = subplot(dim(1),dim(2),i);
    box on; hold on; 
end

%% subfunction : merge subplots
function sax = MergeSubplots(nmrg,dim,ht,i)
    if i == nmrg(1) 
        sax = MakeSubplotDims(dim,nmrg);
    elseif i >= nmrg(2) % if below the big plot, push all other subplots down/over
        sax = MakeSubplotDims(dim,i+1);
    else % make the plot!
        sax = MakeSubplotDims(dim,i);
    end
    if strcmp(dim,'v') % vertically organized plots should have a double-tall merged plot
        ht = ht.*2; 
    end
    ChangeHeight(sax,ht);
end

%% subfunction : make subplot a specified height
function ChangeHeight(ax,ht)
    switch ht
        case 'x' % do nothing
        otherwise 
            ax.Position(end) = ht; % change the height 
    end
end
