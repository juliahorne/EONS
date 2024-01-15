%% ============================= EONS Model =============================== 
% Julia Horne, 2023

% load a series of model runs with the same tuning but different mantle factors
% plot the mantle forcing (constant or evolving from different initial
% multipliers), the oxygen history and the forg history that emerge. 

% An output folder should hold all test runs: load before calling this plot
% function by typing:
% >> load OUTPUT/fman_Tests

function Plot_MantleTest(alltests)
il = 1:1:10;                                     % the values of fman to actually plot

[~,ax] = plot_GenerateSubplots(8,[4 2]);
plot_SqueezeSubplots(ax,'v',2,'tight','n','y');  % smaller subplot spacing
cord = plot_ColorOrder(ax(1),length(il)./6.66,'rainbow');     % make more automatic colors!
v = alltests{1}{1}.v;
for ix = 1:length(alltests)                      % iterate through each test type
   for iy = 1:length(alltests{ix})               % iterate through each fman value
       if any(il == alltests{ix}{iy}.out.inp.fman) % if one wants to skip some values (ie. if il = 1:2:10)
            [t,r,~,flux,~,~,~,~,inp,~,~,~] = UnpackOutput(alltests{ix}{iy}.out);
            tf = TotalOceanFluxes(flux,{'forg','ammon','nitr','mtrophy','burial'}); 
            koxy = tf.burial.OC  ./ ...          % oxygen geologic sources = net burial of reductants
                (0.5.*(flux.volc.CH4 + flux.meta.CH4 + flux.mantle.CH4) ...
                + 0.75.*(flux.volc.NH4 + flux.meta.NH4 + flux.mantle.NH3) ...
                + 0.25.*flux.mantle.FeO + 0.5.*flux.mantle.Fe2SiO4 ); % over sinks = mantle/volcanic/metamorphic reductant outgassing
            % add labels and reference lines only on the first go around 
            if ix == 1 && iy == 1 
                ylabel(ax(1),{'Mantle Reductant';'flux (mol FeO/yr)'});
                lg = legend(ax(1),'-DynamicLegend','location','southwest','NumColumns',5);
                lsz = get(lg,'ItemTokenSize'); lg.ItemTokenSize = lsz./2;
                ylabel(ax(3),'Oxygen level (PAL)'); 
                ylabel(ax(5),{'Organic carbon burial';'fraction (f_{org})'});
                ylabel(ax(7),{'Oxygen sources : sinks'; '(k_{oxy})'})
                for is = 1:length(ax)
                    set(ax(is),'xlim',[1e7 t(end)],'ColorOrder',cord); % set x axis limits and color order for all subplots
                    if is == 1 || is == 2        % only add labels for xlines on top subplots
                        plot_Transitions(ax(is),v,0,'on',10);
                    else
                        plot_Transitions(ax(is),v,0,'off',10);
                    end
                    if is >= length(ax)-1        % remove extranneous x axis labels
                        plot_Agescale(ax(is),'ga');  
                        ax(is).XLabel.FontSize = 10; 
                        ax(is).XTickLabels{end} = '';
                    end
                end
            end
            if iy == 1                            % set specific y limits on first round of plotting
               set(ax(ix),'yscale','lin','ylim',[1e11 3.5e12]);
               set(ax(ix+2),'yscale','log','ylim',[1e-10 1e1]);
               set(ax(ix+4),'yscale','lin','ylim',[0 0.35]);
               set(ax(ix+6),'yscale','log','ylim',[1e-4 20]);
            end
            % now actually plot the output
            if ix == 1                           % only add legend on first plot!
               plot(ax(ix),t,flux.mantle.FeO,'DisplayName',num2str(inp.fman)); % plot mantle influx
            else
               plot(ax(ix),t,flux.mantle.FeO,'HandleVisibility','off');
            end
            % moving down the column 
            plot(ax(ix+2),t,r.a.O2./v.atm.O2pal); % plot oxygen level
            plot(ax(ix+4),t,tf.forg);             % plot f_org
            plot(ax(ix+6),t,koxy); yline(ax(ix+6),1,':'); % plot oxygen balance
            % clear and move on
            clear t tf out r flux koxy
       else
           continue
       end
    end
end

% label all subplots with letters and then print
plot_LabelSubplots(ax,'out','alpha',10);
PrintPDFToFolder(8.25,8.75,'MantleTest',v.figfolder);
end

