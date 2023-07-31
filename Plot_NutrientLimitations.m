%% ============================= EONS Model =============================== 
% Julia Horne, 2023

% Plot the evolution of N and P limitations on productivity and how C:P
% ratio changes in buried organics (proxy for P scavenging efficiency).

function Plot_NutrientLimitations(out,v)

[t,r,conc,flux,~,~,~,~,~,~,~,~] = UnpackOutput(out);
c       = VolumetricConcentrations(r,conc,v);               % calculate concentrations in mol/m3
lim     = Limitations(r,conc,c,v);                          % calculate limitations based on concentrations

%% plot N and P limitations on biosphere
figure(); hold on; box on;
plot(t,lim.P,'color',v.color.H3PO4);
plot(t,lim.N,'color',v.color.fixN);
legend('H_3PO_4','fixed N','location','north','FontSize',10)
set(gca,'FontSize',10,'xlim',[1e8 t(end)],'FontSize',10);
ylabel('Nutrient limitation (L_{assim,x})'); 

% denote bio transitions, change to agescale, and print
plot_Transitions(gca,v,0,'on',10);
plot_Agescale(gca,'ga');
PrintPDFToFolder(7,5.5,'NPLim','/FIGURES');

%% Plot C:P of organic matter burial fluxes (proxy for P scavenging)
figure(); hold on; box on;
plot(t,flux.burial.OC.n./flux.burial.OP.n,'-.'); 
plot(t,flux.burial.OC.z./flux.burial.OP.z,':');
legend('n','z','location','north','fontsize',10)
ylabel('Buried organic matter C:P');
set(gca,'FontSize',10,'xlim',[1e8 t(end)],'ylim',[105 115]);
% add reference line denoting RR
yline(106,':','Modern ratio 106:1','LineWidth',2,'handlevisibility','off',...
    'labelhorizontalalignment','center','FontSize',10);

% denote bio transitions, change to agescale, and print
plot_Transitions(gca,v,0,'on',10);                      
plot_Agescale(gca,'ga');
PrintPDFToFolder(7,5.5,'PScav','/FIGURES');

end