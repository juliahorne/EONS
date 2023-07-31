%% ============================= EONS Model =============================== 
% Julia Horne, 2022

% Inputs are axis handles (ax), standards/parameters (v) containing the
% times for transitions to initiate, time to plot xline separator for
% references (xlinet), labels on/off (lblon) for the lines and font size
% for all labels (fontsz). Everthing after v is optional!

function plot_Transitions(ax,v,xlinet,lblon,fontsz)
    if ~exist('xlinet','var')
       xlinet = 0;
    end
    if ~exist('lblon','var')
       lblon = 'on'; 
    end
    if ~exist('fontsz','var')
       fontsz = 8;  
    end
    % generate labels if each x line is supposed to be labeled
    if strcmp(lblon,'on')
       lbr = 'References'; lls = {'P','B','F','L'};
    else
        lbr = ''; lls = {'','','',''};
    end
    % plot lines with labels for start of photosynthesis, body size growth,
    % fungi, and vascular plants
    times = {'photo','body','fungi','plant'}; 
    for it = 1:length(times)
        PlotXline(ax,v.td.(['init',times{it}]),lls{it},'off','center','top','horizontal',fontsz);
    end
    % delineate where the model runs end and references begin
    if xlinet > 0
        PlotXline(ax,xlinet,lbr,'off','right','middle','aligned',fontsz,'-');
    end
    
end

%% Subfunction: plot an x line with the follwing specifications and labels

% inputs are axis handles for plot or subplot (axs), the value to plot
% (val), the label for the line (label), if the handle should be visible
% in a legend (handlevis = 'on' or 'off'), horizontal and vertical
% alignment and orientation for on-plot label (horz, vert, orient) if 
% desired (if not, enter 'x' for any), label font size and line style 
% (optional entries). 

function xl = PlotXline(axs,val,label,handlevis,horz,vert,orient,fontsize,ls)
   if ~exist('fontsize','var')
        fontsize = 8; 
   end
   if ~exist('ls','var')
       ls = '--';
   end
   xl = xline(axs,val,ls,label,'HandleVisibility',handlevis);             % plot the xline on the axis designated
   xl.LineWidth = 0.5; xl.FontSize = fontsize; 
   switch horz
       case 'x'
       otherwise
           xl.LabelHorizontalAlignment = horz; 
   end
   switch vert
       case 'x'
       otherwise 
           xl.LabelVerticalAlignment = vert;
   end
   switch orient
       case 'x'
       otherwise
           xl.LabelOrientation = orient; 
   end
end 
