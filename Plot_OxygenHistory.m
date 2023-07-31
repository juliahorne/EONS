%% ============================= EONS Model =============================== 
% Julia Horne, 2023

% Comparative plot for model atmospheric oxygen output against oxygen 
% proxies from various literature sources (provided in .csv file). Optional
% input for naming the printed PDF. 

function Plot_OxygenHistory(out,v,pdfname)

if ~exist('pdfname','var')
    pdfname = 'O2ProxyComp';
end

[t,r,~,~,~,~,~,~,~,~,~,~] = UnpackOutput(out);

% Reading from some figures, we estimate duration of different MIF-S
% signals and plot these separately as solid lines (start-finish) and
% arrows (up/down depending on min/max estimates)
addpath ~/EONS/arrow.m   

% import the CSV file of different oxygen level estimates
warning('off','all');                                           % ignore warnings re: column header data IDGAF
oxy = readmatrix('OxygenProxiesTable.csv','Range','A2:D80');    % use readmatrix if only numerical data
txt = readtable('OxygenProxiesTable.csv');                      % use readtable when strings are involved
mya = oxy(:,1);                                                 % first column is ages in Mya
opct= oxy(:,2);                                                 % second column is the "best estimate" oxygen levels at each age in percentages (atm mixing ratio)
minsd= oxy(:,3);                                                % third column is mean-1st std deviation (or general lower boundary estimate)
maxsd= oxy(:,4);                                                % fourth column is mean+1st std deviation (or general upper bounary estimates)

%% organize
for ix = 1:size(mya,1)
    times(ix)  = 4e9 - mya(ix).*1e6;                            % convert mya to time within model run (ie. anything older than 4000 Ma is not plotted)
    oxygen(ix) = opct(ix)./21;                                  % convert atm percentage to PAL
    rngs(ix,:) = [minsd(ix), maxsd(ix)];                        % convert percent std deviations from mean to PAL
    auth	   = char(table2array(txt(ix,5)));                  % take out the author names 
    data{ix}   = char(table2array(txt(ix,6)));                  % take out the type of dataset
    note{ix}   = char(table2array(txt(ix,7)));                  % take out extra notes (max/mins)
    if contains(auth,'Izon')
        symb(ix) = 's'; 
    elseif contains(auth,'Glass')
        symb(ix) = 'o';
    elseif contains(auth,'Farquhar')
        symb(ix) = '*'; 
    elseif contains(auth,'Johnson')
        symb(ix) = 'd'; 
    elseif contains(auth,'Krause')
        symb(ix) = 'x'; 
    elseif contains(auth,'Sperling')
        symb(ix) = '<';
    elseif contains(auth,'Bellefroid')
        symb(ix) = '>';
    elseif contains(auth,'Canfield')                             % == zhang (lead author) who appears as a name in a few of these references
        symb(ix) = '^';
    elseif contains(auth,'Bjerrum')
        symb(ix) = 'h'; 
    elseif contains(auth,'Claire')
        symb(ix) = 'v';
    elseif contains(auth,'Ono')                                  % == Luo (lead author) appears as second author in izon
        symb(ix) = 'p'; 
    end
end

%% now plot!
ax(1) = subplot(3,1,[1 2]); hold on; box on; 
arrowlist = {}; count = 1; 
for ixs = 1:length(times)
    if rngs(ixs,1) == 0
        err(ixs,:) = rngs(ixs,:)./21;                           % convert to PAL
    else
        err(ixs,:) = [oxygen(ixs) - (rngs(ixs,1)./21), (rngs(ixs,2)./21) - oxygen(ixs)]; 
    end
    if strcmp(data{ixs},'MIFS') || contains(data{ixs},'uranite') 
% these proxies provide upper/lower boundaries and the "ranges" are instead start/stop times! Notes column will tell arrow direction
       arrowlist{count} = {[times(ixs), oxygen(ixs), rngs(ixs,:)], note{ixs}};
       count = count + 1; 
       err(ixs,:) = [0 0]; 
    end     
    mksz = changesymbsize(symb(ixs)); 
    eb = errorbar(times(ixs),oxygen(ixs),err(ixs,1),err(ixs,2),...
        'Marker',symb(ixs),'MarkerSize',mksz,'LineWidth',1,...
        'HandleVisibility','off','color',v.color.Archean);
end
% plot model output 
plot(t,r.a.O2./v.atm.O2pal,'color',v.color.O2); 
set(gca,'yscale','log','xscale','lin','xlim',[1e8 t(end)+5e7],'ylim',[1e-10 10]); 
ylabel('Atmospheric O_2 (PAL)'); 

% add arrows to the points that are only max/min estimates
AddMaxMinArrows(arrowlist,10); 

% Final touches!
plot_Transitions(gca,v,0,'on',10);                          % denote when biological transitions occur
plot_Agescale(gca,'ga');                                    % make x axis in Ga
 
%% show a subplot with only the recent, well constrained oxygen history
ax(2) = subplot(3,1,3); hold on; box on;
for ixx = 1:length(times)
    if times(ixx) >= 3.5e9                                  % constrains the phanerozoic (>= 3.5e9 into model run)
        mksz = changesymbsize(symb(ixx)); 
        eb= errorbar(times(ixx),oxygen(ixx),err(ixx,1),err(ixx,2),...
            'Marker',symb(ixx),'MarkerSize',mksz,'LineWidth',1.2,...
            'HandleVisibility','off','color',v.color.Archean);
    end
end
% plot model output again
plot(t,r.a.O2./v.atm.O2pal,'color',v.color.O2);
plot_Transitions(gca,v,0,'on',10);
set(gca,'yscale','lin','xscale','lin','xlim',[3.5e9 t(end)+5e6],'ylim',[0.25 2]); 
ylabel('Atm. O_2 (PAL)'); 
plot_Agescale(gca,'ma');                                    % make x axis in Ma

% add subplot labels (a, b) 
plot_SqueezeSubplots(ax,'v',1,'normal','y','y');            % doesn't squeeze, only widens
plot_LabelSubplots(ax,'out','alpha',10);
PrintPDFToFolder(7,5.5,pdfname,v.figfolder);

end

%% subfunction: add directional arrows 
function AddMaxMinArrows(arwls,arlen)
    for ia = 1:length(arwls)
        updown = arwls{ia}{2};                                  % denotes up/down direction of arrows
        input = arwls{ia}{1};                                   % times, oxygen level, timespan (start, stop)
        xs = 4e9 - ([input(3)  input(4)].*1e6);                 % date around the point
        ys = [input(2)   input(2)];                             % oxygen level 
        plot(xs,ys,'-k'); 
        if strcmp(updown,'min')
           arrow([input(1),input(2)],[input(1),input(2)*arlen],5,50,30); % positions xy, positions xy (higher), arrowhead length, base angle,tip angle
        else % rest are maximums (down arrows)
           arrow([input(1),input(2)],[input(1),input(2)/arlen],5,50,30); % positions xy, positions xy (lower), arrowhead length, base angle,tip angle
        end
    end

end

%% subfunction: change symbol size
function mksz = changesymbsize(sm)
    switch sm
        case {'*','x','h','p'}
            mksz = 10;
        case 's'
            mksz = 12; 
        otherwise
            mksz = 7;
    end

end
