% Reservoir and Flux Comparison Plots
% Julia Horne, 2019
%
% Plots show comparisons of modern-day model run output values for fluxes
% and reservoirs against literature values (present day) as part of the
% model tuning and validation process. Outputs a dual plot of reservoir
% sizes and flux rates, with colored circles/lines denoting the literature
% estimates (compiled in LiteratureComparisonRanges.m) from various
% authors. 


function Plot_output(out,v)

% Compile model output and literature values
LiteratureReference; 
[t,r,~,flux,gasex,~,~,T,~,~,~,~] = UnpackOutput(out);
rf = r; 

box = {'s','d','z','n'}; 
for ib = 1:length(box)
    % calculate fixed N reservoir in ocean boxes
    rf.(box{ib}).fixN = rf.(box{ib}).RN + rf.(box{ib}).HNO3;                   % mol N; fixed N
    rf.(box{ib}).N = 2.*rf.(box{ib}).N2 +  rf.(box{ib}).fixN ;                 % mol N; total marine Nitrogen
end

spe = {'OC','N2','HNO3','RN','DIC','TA','O2','H3PO4','CH4','CaCO3','fixN','FeOH3'}; 
[rf.t] = TotalOceanReservoirs(r,spe);                                        % total ocean reservoirs for these species
% add together reservoirs to create "total ocean" reservoirs
rf.t.LB = rf.s.LB;                                                            % mol C; living marine biosphere
rf.t.POC = rf.s.LB + rf.t.OC ;                                                 % mol C; total organicis (LB + DB = POC) 
rf.t.PON = rf.t.POC / v.const.CNratio;                                        % mol N; N in total marine biomass
% total unreactive + continental (tc)  reservoirs
rf.tc.OC = rf.c.OC + rf.u.OC; 
rf.tc.CaCO3 = rf.c.CaCO3 + rf.u.CaCO3 + rf.o.CaCO3;
rf.tc.NH4 = rf.c.NH4 + rf.o.NH4;
rf.tc.PO4 = rf.c.CP + rf.c.SP + rf.u.SP + rf.u.CP;
rf.tc.Fe = 2.*rf.c.Fe2SiO4 + rf.c.FeOH3 + rf.u.FeOH3 + 2.*rf.c.Fe2O3; 

% combine fluxes into "total ocean" fluxes
fluxes = {'nitr','precip','diss','denitN','denitC','metha','ammon','revweather','burial','axrm','oxrm'};
[flux.t] = TotalOceanFluxes(flux,fluxes); 
flux.fixN = flux.fixation.N2.*2 + flux.fix.newN;                            % mol N; N2 fixation (anoxygenic + oxygenic)

% rename some fluxes for convenience 
flux.prod.totC = flux.prod.CO2;                                             % mol C; total marine productivity    
flux.prod.totN = flux.prod.NH3 + flux.prod.NH4 + flux.prod.HNO3;            % mol N; total N productivity
flux.t.assim = flux.assim ./ v.const.CNratio;                               % mol N; total marine assimilation
flux.t.methox = flux.methox;                                                % mol C; net CH4 oxidation
flux.t.mtrophy = flux.mtrophy; 
flux.t.mantle = flux.mantle.FeO; 
flux.t.Hesc = flux.Hesc.H;
flux.acc.PO4 = flux.acc.SP + flux.acc.CP; 

%% Set up the figure (2x1 plots, reservoirs and fluxes)
[fg,ax] = plot_GenerateSubplots(2,'h'); 
plot_SqueezeSubplots(ax,'h',1,'tight');      

%% Now plot reservoirs
rboxlist = {'t','c','z','d','n','s','a'};
res{1} = {'FeOH3','O2','H3PO4','HNO3','PON','fixN','N2','N','TA','CaCO3','DIC','OC'};  % total ocean
res{2} = {'Fe','PO4','OP','NH4','ON','CaCO3','OC'};                                    % continent + unreactive seds (u)
res{3} = {'FeOH3','O2','H3PO4','HNO3','RN','N2','TA','CaCO3','DIC','OC','CH4'};        % reactive (z+n) sediments
res{4} = {'FeOH3','O2','H3PO4','HNO3','RN','N2','TA','CaCO3','DIC','OC','CH4'};        % deep ocean
res{5} = {'FeOH3','O2','H3PO4','HNO3','RN','N2','TA','CaCO3','DIC','OC','CH4'};        % shallow sediments
res{6} = {'FeOH3','O2','H3PO4','fixN','HNO3','RN','N2','TA','CaCO3','DIC','OC','LB','CH4'};% surace ocean
res{7} = {'NH3','N2','CO2','CH4','O2'};                                     % atmosphere  
Resplotlist = {}; % y axis labels
ik = 1; % counter for y labels

for ibox = 1:length(rboxlist)
    rbox = rboxlist{ibox}; % one reservoir
    spes = res{ibox}; % all of the species fit to plot
    for ispec = 1:length(spes)
        spe = spes{ispec}; % the species to plot
        if strcmp(rbox,'c')
            if isfield(rf.tc,spe) % if a total cont+u sed reservoir doesn't exist, just plot the c reservoir
                plot(ax(1),rf.tc.(spe),ik,'.','MarkerSize',35,'color','k');
            else 
                plot(ax(1),rf.(rbox).(spe),ik,'.','MarkerSize',35,'color','k');
            end % all total cont references are in "c" reservoir
            plot(ax(1),rx.c.(spe),ik,'o','LineWidth',1.5,'MarkerSize',10,'Color',rxcol.c.(spe));
            plot(ax(1),rxr.c.(spe),[ik ik],rsym.c.(spe){1},rsym.c.(spe){2},rsym.c.(spe){3},'color',rrcol.c.(spe),'MarkerSize',7);
        else
            plot(ax(1),rf.(rbox).(spe),ik,'.','MarkerSize',35,'color','k');
            plot(ax(1),rx.(rbox).(spe),ik,'o','LineWidth',1.5,'MarkerSize',10,'Color',rxcol.(rbox).(spe));
            plot(ax(1),rxr.(rbox).(spe),[ik ik],rsym.(rbox).(spe){1},rsym.(rbox).(spe){2},rsym.(rbox).(spe){3},'color',rrcol.(rbox).(spe),'MarkerSize',7);
        end
        Resplotlist{ik} = sprintf([spe,', ',rbox]);
        ik = ik +1;

    end
end
yticks(1:ik-1);
yticklabels(Resplotlist); 
set(gca,'xscale','log','XMinorTick','off','YGrid','on','ylim',[0 ik],'xlim',[1e6 5e22]);
xlabel('Reservoir Size (mol)');


%% Now plot fluxes
fboxlist = {'t','c','u','z','d','n','s'};
fluxes{1} = {'axrm','oxrm','denit','nitr','diss','precip','revweather',...
    'mantle','Hesc','methox'};                                     % total fluxes
fluxes{2} = {'meta','wthr'};                                       % continent
fluxes{3} = {'acc','volc'} ;                                       % unreactive
fluxes{4} = {'mtrophy','ammon','metha','sed','burial'};            % deep sediments
fluxes{5} = {'mtrophy','ammon','metha','export'} ;                 % deep ocean
fluxes{6} = {'mtrophy','ammon','metha','sed','burial'};            % shallow sediments
fluxes{7} = {'mtrophy','ammon','metha','death','ferrotrophy','fix','assim','prod'};% surace ocean

%% 
Fluxplotlist = {}; % y axis labels
ir = 1;  % counter for y labels
for ibox = 1:length(fboxlist)
    rbox = fboxlist{ibox}; % assign box
    if strcmp(rbox,'t')
       tflist = fluxes{ibox};
       for ifx = 1:length(tflist)
           flx = tflist{ifx};
           clist = {v.color.No,v.colorf.C,v.colorf.N,v.colorf.N,v.colorf.C,v.colorf.C,v.colorf.C,v.colorf.Fe,v.colorf.C,v.colorf.O};
           plot(flux.t.(flx)(end),ir,'.','MarkerSize',30,'color','k');
           plot(fluxx.t.(flx),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ifx});
           plotreferencerange(rflux.t.(flx),rsym.t.(flx){1},rsym.t.(flx){2},rsym.t.(flx){3},clist{ifx},ir); 
           Fluxplotlist{ir} = plot_DescriptiveHandle(flx,'t','rirst');  % descriptive legend designation
           ir = ir+1; 
       end 
    end
    flist = fluxes{ibox}; % list fluxes in box 
    for ifx = 1:length(flist)
        flx = flist{ifx} ; % assign flux type
       if strcmp(flx,'prod')
          splst = {'totN','CO2',};%'O2'};
          clist = {v.colorf.N,v.colorf.C};%,v.colorf.O};
          for ispec = 1:length(splst)
              plot(flux.(flx).(splst{ispec})(end),ir,'.','MarkerSize',30,'color','k');
              plot(fluxx.prod.(splst{ispec}),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
              plotreferencerange(rflux.(flx).(splst{ispec}),rsym.(flx).(splst{ispec}){1},rsym.(flx).(splst{ispec}){2},rsym.(flx).(splst{ispec}){3},clist{ispec},ir); 
              Fluxplotlist{ir} = plot_DescriptiveHandle(flx,splst{ispec},'rlast');  % descriptive legend designation
              ir = ir+1; 
          end
         elseif strcmp(flx,'gasex')
             splst = {'NH3','N2','CO2','O2'};
             clist = {v.colorf.N,v.color.No,v.colorf.C};
             sym = {'-^','.','-s'};
             for ispec = 1:length(gaslst)
                  plot(abs(gasex{iy}.(gaslst{ispec})(end)),ir,'.','MarkerSize',35,'color','k');
                  plot(fluxx.(flx).(splst{ispec}),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                  plotreferencerange(rflux.(flx).(splst{ispec}),rsym.(flx).(splst{ispec}){1},rsym.(flx).(splst{ispec}){2},rsym.(flx).(splst{ispec}){3},clist{ispec},ir); 
                  Fluxplotlist{ir} = plot_DescriptiveHandle(flx,gaslst{ispec},'rlast');  % descriptive legend designation
                  ir = ir+1; 
             end
       elseif strcmp(flx,'death')
              plot(flux.(flx).OC(end),ir,'.','MarkerSize',35,'color','k');
              plot(fluxx.(flx)(end),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',v.colorf.C);
              plotreferencerange(rflux.(flx),rsym.(flx){1},rsym.(flx){2},rsym.(flx){3},v.colorf.C,ir); 
              Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{},'rlast');               % descriptive legend designation
              ir = ir+1; 
       elseif strcmp(flx,'fix')
              plot(flux.fixN(end),ir,'.','MarkerSize',35,'color','k');
              plot(fluxx.(flx)(end),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',v.colorf.N);
              plotreferencerange(rflux.(flx),rsym.(flx){1},rsym.(flx){2},rsym.(flx){3},v.colorf.N,ir); 
              Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{},'rlast');  % descriptive legend designation
              ir = ir+1; 
       elseif strcmp(flx,'assim')
              plot(flux.t.(flx)(end),ir,'.','MarkerSize',35,'color','k');
              plot(fluxx.(flx)(end),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',v.colorf.N);
              plotreferencerange(rflux.(flx),rsym.(flx){1},rsym.(flx){2},rsym.(flx){3},v.colorf.N,ir); 
              Fluxplotlist{ir} = plot_DescriptiveHandle(flx,'FixedN','rlast');  % descriptive legend designation
              ir = ir+1; 
       elseif strcmp(flx,'ferrotrophy')
              plot(flux.(flx).CO2(end),ir,'.','MarkerSize',35,'color','k');
              plot(fluxx.(flx).CO2(end),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',v.colorf.C);
              plotreferencerange(rflux.(flx).CO2,rsym.(flx).CO2{1},rsym.(flx).CO2{2},rsym.(flx).CO2{3},v.colorf.C,ir); 
              Fluxplotlist{ir} = plot_DescriptiveHandle(flx,'CO2','rlast');  % descriptive legend designation
              ir = ir+1; 
      elseif strcmp(flx,'ammon')
            splst = {'H3PO4','RN','OC'};
            if strcmp(rbox,'s') || strcmp(rbox,'d')
                clist = {v.color.No,v.color.No,v.colorf.C};
                rlist = {v.color.No,v.colorf.N,v.colorf.C};
            elseif strcmp(rbox,'n') 
                clist = {v.color.No,v.color.No,v.colorf.C};
                rlist = {v.color.No,v.color.No,v.colorf.C};
            else
                clist = {v.colorf.P,v.color.No,v.color.No};
                rlist = {v.colorf.P,v.color.No,v.color.No};
            end
            sym = {'.','-<','->'};
            for ispec = 1:length(splst)
                plot(flux.(flx).(splst{ispec}).(rbox)(end),ir,'.','MarkerSize',35,'color','k');
                plot(fluxx.(flx).(splst{ispec}).(rbox),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}).(rbox),rsym.(flx).(splst{ispec}).(rbox){1},rsym.(flx).(splst{ispec}).(rbox){2},rsym.(flx).(splst{ispec}).(rbox){3},rlist{ispec},ir); 
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{rbox,splst{ispec}},'rlast');  % descriptive legend designation

                ir = ir+1; 
            end 
      elseif strcmp(flx,'metha')
            splst = {'H3PO4','RN','OC'};
            clist = {v.color.No,v.color.No,v.color.No};
            for ispec = 1:length(splst)
                plot(flux.(flx).(splst{ispec}).(rbox)(end),ir,'.','MarkerSize',35,'color','k');
                plot(fluxx.(flx).(splst{ispec}).(rbox),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}).(rbox),rsym.(flx).(splst{ispec}).(rbox){1},rsym.(flx).(splst{ispec}).(rbox){2},rsym.(flx).(splst{ispec}).(rbox){3},clist{ispec},ir); 
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{rbox,splst{ispec}},'rlast');  % descriptive legend designation
                ir = ir+1; 
            end
      elseif strcmp(flx,'mtrophy')
            plot(flux.(flx).(rbox)(end),ir,'.','MarkerSize',35,'color','k');
            plot(fluxx.(flx).(rbox),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',v.color.No);
            plotreferencerange(rflux.(flx).(rbox),rsym.(flx).(rbox){1},rsym.(flx).(rbox){2},rsym.(flx).(rbox){3},v.color.No,ir); 
            Fluxplotlist{ir} = plot_DescriptiveHandle(flx,rbox,'rlast');                 % descriptive legend designation
            ir = ir+1; 
         elseif strcmp(flx,'wthr')
            splst = {'carb','sil','oxi','PO4','NH4','OP','ON'};
            clist = {v.colorf.C,v.colorf.C,v.colorf.C,v.colorf.P,v.colorf.N,v.color.No,v.colorf.N};
            for ispec = 1:length(splst)
                plotflux(ax(2),flux,fluxx,rflux,rsym,flx,splst{ispec},ir,clist{ispec}); % do plotting
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,splst{ispec},'rlast');  % descriptive legend designation
                ir = ir+1; 
            end 
        elseif strcmp(flx,'volc')
            splst = {'NH4','OC','carb'};
            clist = {v.colorf.N, v.color.No, v.colorf.C}; 
            for ispec = 1:length(splst)
                plot(flux.volc.(splst{ispec})(end),ir,'.','MarkerSize',35,'color','k');
                plot(fluxx.(flx).(splst{ispec}),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}),rsym.(flx).(splst{ispec}){1},rsym.(flx).(splst{ispec}){2},rsym.(flx).(splst{ispec}){3},clist{ispec},ir);
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,splst{ispec},'rlast');  % descriptive legend designation

                ir = ir+1; 
            end
        elseif strcmp(flx,'meta')
            splst = {'NH4','OC','CaCO3'};
            clist = {v.color.No, v.colorf.C, v.colorf.C}; 
            for ispec = 1:length(splst)
                plot(flux.(flx).(splst{ispec})(end),ir,'.','MarkerSize',35,'color','k');
                plot(fluxx.(flx).(splst{ispec}),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}),rsym.(flx).(splst{ispec}){1},rsym.(flx).(splst{ispec}){2},rsym.(flx).(splst{ispec}){3},clist{ispec},ir);
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,splst{ispec},'rlast');  % descriptive legend designation

                ir = ir+1; 
           end

         elseif strcmp(flx,'export')
             splst = {'CaCO3','OC'};
             clist = {v.colorf.C,v.colorf.C};
             for ispec = 1:length(splst)
                 for iy = 1:length(out.r)
                    plot(flux.(flx).(splst{ispec})(end),ir,'.','MarkerSize',35,'color','k');
                 end
                plot(fluxx.(flx).(splst{ispec}),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}),rsym.(flx).(splst{ispec}){1},rsym.(flx).(splst{ispec}){2},rsym.(flx).(splst{ispec}){3},clist{ispec},ir);
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,splst{ispec},'rlast');  % descriptive legend designation
                ir = ir+1; 
             end
     
       elseif strcmp(flx,'sed')
            splst = {'CaCO3','OC'};
            clist = {v.colorf.C,v.colorf.C};
            for ispec = 1:length(splst)
                for iy = 1:length(out.r)
                    plot(flux.(flx).(splst{ispec}).(rbox)(end),ir,'.','MarkerSize',35,'color','k');
                end
                plot(fluxx.(flx).(splst{ispec}).(rbox),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}).(rbox),rsym.(flx).(splst{ispec}).(rbox){1},rsym.(flx).(splst{ispec}).(rbox){2},rsym.(flx).(splst{ispec}).(rbox){3},clist{ispec},ir);
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{rbox,splst{ispec}},'rlast');  % descriptive legend designation

                ir = ir+1; 
            end
            
        elseif strcmp(flx,'diff')
            splst = {'H3PO4','HNO3','RN','N2','DIC','O2'};
            flist = {v.color.No,v.color.No,v.color.No,v.color.No,v.colorf.C,v.color.No};
            for ispec = 1:length(splst)
                for iy = 1:length(out.r)
                    plot(abs(flux.(flx).(splst{ispec}).(rbox)(end)),ir,'.','MarkerSize',35,'color','k');
                end
                plotreferencerange(rflux.(flx).(splst{ispec}).(rbox),rsym.(flx).(splst{ispec}).(rbox){1},rsym.(flx).(splst{ispec}).(rbox){2},rsym.(flx).(splst{ispec}).(rbox){3},clist{ispec},ir);
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{rbox,splst{ispec}},'rlast');  % descriptive legend designation

                ir = ir+1; 
            end
        elseif strcmp(flx,'burial')
            splst = {'FeOH3','PO4','OP','ON','CaCO3','OC'};
            clist = {v.colorf.Fe,v.colorf.P,v.color.No,v.colorf.N,v.colorf.C,v.colorf.C};
            rlist = {v.colorf.Fe,v.colorf.P,v.colorf.P,v.colorf.N,v.colorf.C,v.colorf.C};
            for ispec = 1:length(splst)
                for iy = 1:length(out.r)
                    plot(flux.(flx).(splst{ispec}).(rbox)(end),ir,'.','MarkerSize',35,'color','k');
                end
                plot(fluxx.(flx).(splst{ispec}).(rbox),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',clist{ispec});
                plotreferencerange(rflux.(flx).(splst{ispec}).(rbox),rsym.(flx).(splst{ispec}).(rbox){1},rsym.(flx).(splst{ispec}).(rbox){2},rsym.(flx).(splst{ispec}).(rbox){3},rlist{ispec},ir);
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,{rbox,splst{ispec}},'rlast');  % descriptive legend designation
                ir = ir+1; 
            end
        elseif strcmp(flx,'acc')
            splst = {'PO4','OP','ON','carb','OC'};
            for ispec = 1:length(splst)
                for iy = 1:length(out.r)
                    plot(flux.(flx).(splst{ispec})(end),ir,'.','MarkerSize',35,'color','k');
                end
                Fluxplotlist{ir} = plot_DescriptiveHandle(flx,splst{ispec},'rlast');  % descriptive legend designation
                ir = ir+1; 
            end
       else
           continue
        end
    end
end

yticks(1:ir-1);
yticklabels(Fluxplotlist); xlabel('Fluxes (mol/yr)'); 
set(gca,'xscale','log','XMinorTick','off','YGrid','on','ylim',[0 ir],'xlim',[5e3 1e17]);

%% Final touches
st= sgtitle({['Total runtime: ', num2str(t(end),'%1.e'),' yrs'];...
    ['T_{sur} = ', num2str(T(end),'%.f'), ' K']}); 
st.FontSize = 12;
set(fg,'Units','inches','PaperUnits','inches','PaperSize',pagesz);
PrintPDFToFolder(pagesz(1),pagesz(2),'FinalOutput','/EONSFigs');

end

%% Subfunction : plot flux and refereces
function plotflux(ax,flux,fluxx,rflux,rsym,FNM,SPNM,ir,COLOR) 
switch SPNM
    case 'x'
        plot(ax,flux.(FNM)(end),ir,'.','MarkerSize',35,'color','k');
        plot(ax,fluxx.(FNM),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',COLOR);
        plotreferencerange(ax,rflux.(FNM),rsym.(FNM){1},rsym.(FNM){2},rsym.(FNM){3},COLOR,ir); 
    otherwise % two layers deep
        plot(ax,flux.(FNM).(SPNM)(end),ir,'.','MarkerSize',35,'color','k');
        plot(ax,fluxx.(FNM).(SPNM),ir,'o','LineWidth',1.5,'MarkerSize',10,'Color',COLOR);
        plotreferencerange(ax,rflux.(FNM).(SPNM),rsym.(FNM).(SPNM){1},rsym.(FNM).(SPNM){2},rsym.(FNM).(SPNM){3},COLOR,ir); 
end
end

%% subfunction to plot reference range for each flux
function plotreferencerange(FLUX,SYM,MARKERFACE,FACECOLOR,specolor,ir)
% use the inputted plot commands to plot the reference range
if strcmp(FACECOLOR,'color') % this designation tells me to just use the species color
    plot(ax,FLUX,[ir ir],SYM,'color',specolor,MARKERFACE,specolor,'MarkerSize',7); 
else % use the color designated (helps identify referenced authors!)
    plot(ax,FLUX,[ir ir],SYM,'color',specolor,MARKERFACE,FACECOLOR,'MarkerSize',7);
end

end

