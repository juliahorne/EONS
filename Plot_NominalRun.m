%% ============================= EONS Model =============================== 
% Julia Horne, 2022

% Produces the figures from EONS paper (TBD) results section- reservoirs,
% fluxes of a particular system through time. 
% Optional input for what figure to produce (figs) - note the options below
% in the 'figfocus' switch command. No input for figs will produce 'paper'.

function Plot_NominalRun(out,v,figs)
    if ~exist('figs','var')
        figs = 'paper';
    end
    xlimt = 4.4e9;     % extra length of x axis to allow inclusion of references
    switch figs
        case 'paper' % all the figures in the published paper (with a few extras)
            figfocus = {'geocarb','geoN','geoP','geoOC','bio','ocean','osys','temprf','atm','forcings','mantlesurf'};
        case 'all'
            figfocus = {'balance','bseN','subzone','Silcyc','geocarb','pH','geoN','geoP','geoOC','bio','bioN','ocean','osys','temprf','atm','forcings'};
        case 'nonSS' % show an extra billion years that is definitely NOT at steady state
            figfocus = {'geocarb','geoN','geoP','geoOC','bio','ocean','osys','temprf','atm','forcings'};
            xlimt = 5e9; 
        otherwise
            figfocus = {figs}; 
    end
    xlinet = 4e9; reft = 4.1e9; % when to plot reference x line and start of references
    LiteratureReference; 
	[t,r,conc,flux,gasex,mix,rf,T,inp,tdep,~,~] = UnpackOutput(out);
    cx = VolumetricConcentrations(r,conc,v); 
    % calculate total ocean reservoirs and total ocean fluxes
    spe = {'OC','ON','OP','N2','HNO3','RN','DIC','TA','O2','H3PO4','CH4','CaCO3','fixedN','FeOH3'}; 
    [tr] = TotalOceanReservoirs(r,spe); 
    fxs = {'nitr','precip','diss','denitN','denitC','metha','ammon','revweather','burial','forg','mtrophy','bureff'};
    [tf] = TotalOceanFluxes(flux,fxs); 
    bx = {'s','d','n','z'}; 
    for ib = 1:length(bx)
        r.(bx{ib}).fixN = r.(bx{ib}).RN + r.(bx{ib}).HNO3; 
        r.(bx{ib}).NP   = r.(bx{ib}).fixN ./ r.(bx{ib}).H3PO4; 
    end
    [totalres,~,~,~,~] = SumAllSpecies(r,conc,inp.indx,inp,v);  % calculate total elemental masses
    
    %% Use AGU conventions for font size, figure dimensions, etc. 
    fontsz = 9;                 % text should be >= pt 8 font
    qtrpg  = [3.74 4.52];       % in; quarter-page figure dimensions (width, height; converted from mm)
    fullpg = [7.48 9.055];      % in; full-page figure dimensions 
    % subplot height is relative, should be stretched to fit full and
    % quarter pages.
    spht = 0.13;% this works best because of the squeeze subplots extension!
    % plot references at the end of the line, just after the end of model
    % run but to the left of the outline of the plot
    refx = reft - 0.04e9; %  start around 4.01 or 5.01 billion years instead  of 4, 5 in the case of multiple references
   
    % placeholder references and reference ranges, for unconstrained values
    pref = 1e-15; prng = [1 1].*1e-15; 
    
    %% get to the good plotting shit
    % make printing names and define transition timeframes we will highlight!
    frames = [1e8 xlimt];
    % designate linestyles for ocean boxes as different from default
    % (s,n,d,z) and for cont/seds (u,o,c)
    lsty.s = '-'; lsty.d = '--'; lsty.n = '-.'; lsty.z = ':'; 
    lsty.u = ':'; lsty.o = '--'; lsty.c = '-'; 
    % generate 4 transect plots, and three frames for each selected output
    % (ie. 12 plots!)
    for ig = 1:length(figfocus)
        clear xv c yl nm lbound ubound dim nmrg lSty
        pagesz = fullpg; 
        switch figfocus{ig}
            case 'atm'     % atmospheric species in mixing ratios on ONE axis
                xv = {'O2','N2','CO2','CH4','NH3'}; 
                fg = figure(); ax = axes(); hold on; box on; 
                for ix = 1:length(xv)
                    dispnm = plot_DescriptiveHandle('SPECIES',xv{ix});   % print out a nice stoichiometric name for the species
                    plot(ax,t,r.a.(xv{ix})./v.atm.mol,'color',v.color.(xv{ix}),...
                        'DisplayName',dispnm);
                    PlotRefRange(ax,reft,v.atm.([xv{ix},'pal'])./v.atm.mol,prng,v.color.(xv{ix})); % plot reference and range
                end
                set(gca,'FontSize',10,'xscale','lin','xlim',[1e8 4.25e9],'yscale','log');
                plot_Transitions(gca,v,xlinet,'on',fontsz);
                ylims = get(gca,'ylim'); set(gca,'ylim',[1e-12 1e2]); 
                legend('-DynamicLegend','location','north','NumColumns',3);
                plot_Agescale(ax,'ga');  
                ylabel('Atmospheric relative mixing ratio (mol/mol_{atm})');

            otherwise  % plot on subplots!
            switch figfocus{ig}
                case 'ocean'
                    lgpos = 'southwest'; dim = [5 2];  
                    xx = {'O2','CH4','RN','HNO3','H3PO4','OC','DIC','TA','pH','omega'}; 
                    lbound = [pref 1e-10 1e-5 1e-5 1e-5 1e-8 1e0 1e0 7 0]; ubound = [1e1 1e-2 1e-1 1e-1 1e-2 1e4 2e1 2e1 8.25 3.5]; 
                    for iv = 1:length(xx)
                        dispnm = plot_DescriptiveHandle('SPECIES',xx{iv});
                        bs = {'s','n','d','z'};
                        for ib = 1:length(bs)
                            lSty{iv}{ib} = lsty.(bs{ib});               % assign ocean box line style
                            switch xx{iv}
                                case {'pH','omega'}
                                    switch xx{iv}
                                        case 'omega'
                                             xv{iv}{ib} = conc.Oc.(bs{ib}); c{iv}{ib}  = v.color.CaCO3;
                                             ref{iv}{ib}= pref; rng{iv}{ib} = prng; 
                                        otherwise
                                             xv{iv}{ib} = conc.(bs{ib}).(xx{iv}); c{iv}{ib}  = v.color.FeO;
                                             if ib == 1
                                                 ref{iv}{ib}= 8.1; rng{iv}{ib} = [7.8 8.2]; 
                                             else
                                                 ref{iv}{ib}= pref; rng{iv}{ib} = prng; 
                                             end
                                    end
                                    nm{iv}{ib} = bs{ib}; 
                                    yl{iv} = dispnm; 
                                otherwise
                                    xv{iv}{ib} = cx.(bs{ib}).(xx{iv}); 
                                    nm{iv}{ib} = bs{ib}; 
                                    c{iv}{ib}  = v.color.(xx{iv});
                                    yl{iv} = {['[',dispnm,'] (mol/m^3)']};                                        
                                    if rx.(bs{ib}).(xx{iv}) ~= pref
                                        ref{iv}{ib} = rx.(bs{ib}).(xx{iv})./v.oc.vol.(bs{ib});
                                    else
                                        ref{iv}{ib} = rx.(bs{ib}).(xx{iv}); 
                                    end
                                    if rxr.(bs{ib}).(xx{iv})(1) ~= pref
                                        rng{iv}{ib} = rxr.(bs{ib}).(xx{iv})./v.oc.vol.(bs{ib});
                                    else 
                                        rng{iv}{ib} = rxr.(bs{ib}).(xx{iv});
                                    end
                            end
                        end
                    end
                case 'pH'
                    lgpos = 'south'; dim = [2 2]; 
                    bs = {'s','n','d','z'}; 
                    for ib = 1:length(bs)
                       xv{1}{ib} = conc.(bs{ib}).pH;    c{1}{ib} = v.color.pH;   nm{1}{ib} = bs{ib}; lbound(1) = 6.8; ubound(1) = 8.5;
                       xv{2}{ib} = flux.precip.(bs{ib});c{2}{ib} = v.color.CO3;  nm{2}{ib} = bs{ib};lbound(2) = 1e11;ubound(2) = 5e14;
                       xv{3}{ib} = conc.Oc.(bs{ib});    c{3}{ib} = v.color.CO3;  nm{3}{ib} = bs{ib}; lbound(3) = 0;   ubound(3) = 2;
                       xv{4}{ib} = flux.diss.(bs{ib});  c{4}{ib} = v.color.CO3;  nm{4}{ib} = bs{ib}; lbound(4) = 1e8;ubound(4) = 5e14;
                       for il = 1:length(xv)
                           lSty{il}{ib} = lsty.(bs{ib}); 
                       end
                       switch bs{ib} 
                           case 's'
                               ref{1}{ib} = 8.1;        ref{2}{ib} = pref;           ref{3}{ib} = pref;  ref{4}{ib} = pref;
                               rng{1}{ib} = [7.8 8.2];  rng{2}{ib} = rflux.t.precip; rng{3}{ib} = prng;  rng{4}{ib} = rflux.t.diss; 
                           otherwise
                               for ii = 1:4
                                   ref{ii}{ib} = pref; rng{ii}{ib} = prng;
                               end
                       end
                    end
                    yl{1} = {'Ocean/sediment pH'}; yl{3} = {'Saturation (\Omega_{cal})'}; 
                    yl{2} = {'Precipitation';'(mol C/yr)'}; yl{4} = {'Dissolution';'(mol C/yr)'}; 

                case 'subzone' % c reservoirs, fluxes; N r,f; P r, f; mantle res, fluxes;
                    lgpos = 'south'; dim = [4 2];
                    lbound(1) = 1e15; ubound(1) = 5e22;                                 lbound(3) = 1e14; ubound(3) = 5e20;
                    xv{1} = {r.c.CaCO3, r.c.OC, r.u.OC, r.o.CaCO3+r.u.CaCO3};           xv{3} = {r.c.NH4, r.c.ON, r.u.ON, r.o.NH4}; 
                    ref{1}= {rx.c.CaCO3, rx.c.OC, pref, pref};                          ref{3}= {rx.c.NH4, rx.c.ON, rx.u.ON, rx.o.NH4}; 
                    rng{1}= {rxr.c.CaCO3,rxr.c.OC,prng,prng};                           rng{3}= {rxr.c.NH4,rxr.c.ON, rxr.u.ON, rxr.o.NH4}; 
                    c{1}  = {v.color.CaCO3, v.color.OC, v.color.C, v.color.CO2};        c{3}  = {v.color.NH4, v.color.ON, v.color.N, v.color.N2}; 
                    nm{1} = {'cCaCO_3','cOrg C','uOrg C','u+oCaCO_3'};                  nm{3} = {'cNH_4','cOrg N','uOrg N','oNH_4'};
                    yl{1} = {'Crustal carbon';'reservoirs (mol C)'};                    yl{3} = {'Crustal nitrogen';'reservoirs (mol N)'}; 

                    lbound(5) = 1e14; ubound(5) = 5e20;                                 lbound(2) = 1e10; ubound(2) = 5e13;
                    xv{5} = {r.c.SP, r.c.CP, r.c.OP};                                   xv{2} = {flux.subduct.CaCO3+flux.subduct.carb+flux.subduct.OC,flux.acc.CaCO3+flux.acc.carb+flux.acc.OC,flux.volc.CaCO3+flux.volc.carb+flux.volc.OC}; 
                    ref{5}= {rx.c.SP,rx.c.CP, rx.c.OP};                                 ref{2}= {pref,pref,fluxx.volc.carb}; 
                    rng{5}= {rxr.c.SP,rxr.c.CP, rxr.c.OP};                              rng{2}= {prng,prng,rflux.volc.carb}; 
                    c{5}  = {v.color.SP, v.color.CP, v.color.OP};                       c{2}  = {v.color.C,v.color.CaCO3,v.color.CO2}; % use same colors for all the below subplots
                    nm{5} = {'silicate','carbonate','organic'};                         nm{2} = {'subduction','accretion','volcanism'};
                    yl{5} = {'Crustal phosphate';'reservoirs (mol P)'};                 yl{2} = {'Carbon';'fluxes (mol C/yr)'};

                    lbound(4) = 1e9; ubound(4) = 1e12;                                 lbound(6) = 1e8; ubound(6) = 5e12;
                    xv{4} = {flux.subduct.ON+flux.subduct.NH4,flux.acc.ON+flux.acc.NH4,flux.volc.ON+flux.volc.NH4,flux.cryst.NH4};        xv{6} = {flux.subduct.OP+flux.subduct.CP+flux.subduct.SP,flux.acc.SP+flux.acc.CP+flux.acc.OP,flux.cryst.CP+flux.cryst.SP+flux.cryst.OP}; 
                    ref{4}= {fluxx.subduct.NH4+fluxx.subduct.ON,fluxx.acc.NH4+fluxx.acc.ON,fluxx.volc.NH4+fluxx.volc.ON,fluxx.cryst.NH4}; ref{6}= {pref,pref,pref}; 
                    rng{4}= {rflux.subduct.NH4+rflux.subduct.ON,rflux.acc.NH4+rflux.acc.ON,rflux.volc.NH4+rflux.volc.ON,rflux.cryst.NH4}; rng{6}= {prng,prng,prng}; 
                    c{4}  = {v.color.NH4,v.color.ON,v.color.N,v.color.N2};             c{6} = {v.color.OP,v.color.CP,v.color.SP}; 
                    nm{4} = {'subduction','accretion','volcanism','crystallization'};  nm{6} = {'subduction','accretion','crystallization'};
                    yl{4} = {'Nitrogen';'fluxes (mol N/yr)'};                          yl{6} = {'Phosphorus';'fluxes (mol P/yr)'};

                    lbound(8) = 1e7; ubound(8) = 5e12;                                 lbound(7) = 1e20; ubound(7) = 5e22;           
                    xv{8} = {flux.mantle.Fe2SiO4, flux.mantle.C, flux.mantle.N, flux.mantle.P,flux.mantle.SiO3}; xv{7} = {r.m.N, r.m.C, r.m.P};
                    ref{8}= {pref,pref,pref,pref,pref};                                 ref{7}= {v.m.N,v.m.C,pref}; 
                    rng{8}= {prng,prng,prng,prng,prng};rng{7}= {[3 10].*v.atm.PAN,prng,prng}; 
                    c{8}  = {v.color.Fe,v.color.C,v.color.N,v.color.P,v.color.Archean}; c{7}  = {v.color.N, v.color.C, v.color.P};
                    nm{8} = {'Fe_2SiO_4','CO_2/CH_4','N_2/NH_3','SP','SiO_3',};         nm{7} = {'N','C','P'}; 
                    yl{8} = {'Mantle flux';'(mol/yr)'};                                 yl{7} = {'Mantle reservoir';'(mol)'};

                case 'bio'
                    lgpos = 'south'; dim = 'v'; nmrg = [1 2];
                    lbound(2) = 1e5; ubound(2) = 5e16;                                  lbound(1) = 1e5; ubound(1) = 1e16;
                    xv{2} = {flux.assim, flux.ferrotrophy.CO2};                         xv{1} = {r.s.RN, r.s.HNO3, r.s.RN+r.s.HNO3,r.s.H3PO4,r.s.OC,r.s.FeO};
                    ref{2}= {fluxx.prod.CO2,pref};                                      ref{1}= {rx.s.RN, rx.s.HNO3,rx.s.fixN,rx.s.H3PO4,rx.s.OC,pref}; 
                    rng{2}= {rflux.prod.CO2,prng};                                      rng{1}= {rxr.s.RN,rxr.s.HNO3,rxr.s.fixN,rxr.s.H3PO4,rxr.s.OC,prng};
                    c{2}  = {v.color.CO2, v.color.FeO};                                 c{1}  = {v.color.RN,v.color.HNO3,v.color.fixN,v.color.H3PO4,v.color.OC,v.color.FeO};                                             
                    nm{2} = {'photosynthesis','photo-ferrotrophy'};                     nm{1} = {'NH_3 + NH_4^+','HNO_3','fixed N','H_3PO_4','POC','FeO'}; 
                    yl{2} = {'Primary';'Production (mol C/yr)'};                        yl{1} = {'Surface ocean reservoir (mol)'}; 

                    lbound(3) = 1e5; ubound(3) = 1e14;                                  lbound(4) = 1e10; ubound(4) = 5e16;                                                                             
                    xv{3} = {flux.fixation.N2./2, flux.fix.newN};                       xv{4} = {tf.ammon, tf.denitC, tf.metha};                                                       
                    ref{3}= {pref,fluxx.fix};                                           ref{4}= {fluxx.t.oxrm,fluxx.t.denit.*(106./84.8),pref};         
                    rng{3}= {prng,rflux.fix};                                           rng{4}= {rflux.t.oxrm,rflux.t.denit.*(106./84.8),prng};                                           
                    c{3}  = {v.color.N2,v.color.NH3};                                   c{4}  = {v.color.CO2,v.color.HNO3,v.color.CH4};                                                
                    nm{3} = {'anoxygenic','oxygenic'};                                  nm{4} = {'ammonification','denitrification','methanogenesis'};                                  
                    yl{3} = {'N_2 Fixation';'(mol N/yr)'};                              yl{4} = {'Remineralization';'(mol C/yr)'};                          

                    lbound(5) = 0; ubound(5) = 0.5;                                     lbound(6) = 1e-2; ubound(6) = 1;
                    xv{5} = {tf.forg};                                                  xv{6} = {100.*tf.bureff}; 
                    ref{5}= {0.2};                                                      ref{6}= {0.2}; % Crockford et al. 2023 estimate for modern burial  efficiency percent
                    rng{5}= {prng};                                                     rng{6}= {prng};
                    c{5}  = {v.color.OC};                                               c{6}  = {v.color.OC}; 
                    nm{5} = {''};                                                       nm{6} = {''};
                    yl{5} = {'Organic carbon';'burial fraction (f_{org})'};             yl{6} = {'Organic carbon';'burial efficiency (%)'};

                case 'temprf'
                    lgpos = 'west'; dim = 'v'; nmrg = [1 2]; 
                    lbound(2) = 1e-2; ubound(2) = 100;                                  lbound(1) = 1e-12; ubound(1) = 1e1;
                    xv{2} = {rf.CO2, rf.NH3, rf.CH4};                                   xv{1} = {r.a.CO2./v.atm.mol, r.a.NH3./v.atm.mol, r.a.CH4./v.atm.mol, r.a.N2./v.atm.mol, r.a.O2./v.atm.mol};    
                    ref{2}= {35, 0, 2};                                                 ref{1}= {v.atm.CO2pal/v.atm.mol, v.atm.NH3pal/v.atm.mol, v.atm.CH4pal/v.atm.mol, v.atm.N2pal/v.atm.mol, v.atm.O2pal/v.atm.mol}; 
                    rng{2}= {[32 37],prng,prng};                                        rng{1}= {prng,prng,prng,prng,prng}; 
                    c{2}  = {v.color.CO2, v.color.NH3, v.color.CH4};                    c{1}  = {v.color.CO2,v.color.NH3,v.color.CH4,v.color.N2,v.color.O2};
                    nm{2} = {'CO_2','NH_3','CH_4'};                                     nm{1} = {'CO_2','NH_3','CH_4','N_2','O_2'};         
                    yl{2} = {'Radiative forcing';'(W/m^2)'};                            yl{1} = {'Relative mixing ratio (mol/mol_{atm})'}; 

                    lbound(3) = 1e2; ubound(3) = 310;                                   lbound(4) = 1e3;  ubound(4) = 1500;                                            
                    xv{3} = T;                                                          xv{4} = rf.sc;     
                    ref{3}= 289;                                                        ref{4}= v.S.Pref; 
                    rng{3}= [287 293];                                                  rng{4}= [1360 1362]; 
                    c{3}  = v.color.T;                                                  c{4}  = v.color.Modern; 
                    nm{3} = '';                                                         nm{4} = '';
                    yl{3} = {'Surface';'Temperature (K)'};                              yl{4} = {'Solar Constant';'(W/m^2)'};   

                case 'bseN' % references from Johnson+ Goldblatt 2015 table 13
                    lgpos = 'west'; dim = 'v'; lbound = [1e-5 0 0 0]; ubound(1:4) = 10; 
                    xv{1} = {(r.u.ON+r.c.ON)./v.atm.PAN,(r.o.NH4+r.c.NH4)./v.atm.PAN,r.m.N./v.atm.PAN};xv{2} = {r.u.ON./v.atm.PAN, r.c.ON./v.atm.PAN};  
                    ref{1}= {0.41+1.55, 0.21+0.55, 7.2};                                ref{2}= {0.41, 1.55};
                    rng{1}= {[1.14 2.78],[0.475 1.045], [1.3 13.1]};                    rng{2}= {[0.21 0.61], [0.93 2.17]};
                    c{1}  = {v.color.ON, v.color.NH4, v.color.N};                       c{2}  = {v.color.ON, v.color.N};
                    nm{1} = {'ON','NH_4^+','N'};                                        nm{2} = {'unreactive sed','continent'}; 
                    yl{1} = {'BSE reservoirs';'(PAN)'};                                 yl{2} = {'Organic N';'reservoirs (PAN)'}; 

                    xv{3} = {r.o.NH4./v.atm.PAN, r.c.NH4./v.atm.PAN};                   xv{4} = {(2.*r.a.N2)./v.atm.PAN};  
                    ref{3}= {0.21, 0.55};                                               ref{4}= {1};
                    rng{3}= {[0.195 0.225],[0.28 0.82]};                                rng{4}= {prng};
                    c{3}  = {v.color.NH4, v.color.N};                                   c{4}  = {v.color.N2};
                    nm{3} = {'ocean crust','continent'};                                nm{4} = {''}; 
                    yl{3} = {'NH_4^+ reservoirs';'(PAN)'};                              yl{4} = {'Atmospheric N_2';'(PAN)'}; 

                case 'bioN'
                    lgpos = 'north'; ubound = [5e17 1e16];  dim = 'v';
                    lbound(1) = 1e13;                                                   lbound(2) = 1e7; 
                    xv{1} = {tr.RN, tr.HNO3, tr.fixedN};                                xv{2} = {flux.fix.newN, flux.fixation.N2.*2, tf.nitr, tf.denit}; 
                    ref{1}= {rx.t.RN, rx.t.HNO3, rx.t.fixN};                            ref{2}= {fluxx.fix, pref, fluxx.t.nitr, fluxx.t.denit}; 
                    rng{1}= {rxr.t.RN, rxr.t.HNO3, rxr.t.fixN};                         rng{2}= {rflux.fix, prng, rflux.t.nitr, rflux.t.denit}; 
                    c{1}  = {v.color.RN, v.color.HNO3, v.color.fixN};                   c{2}  = {v.color.NH3, v.color.N2, v.color.NH4, v.color.HNO3};
                    nm{1} = {'NH_3 + NH_4^+','HNO_3','fixed N'};                        nm{2} = {'oxygenic fixation','anoxygenic fixation','nitrification','denitrification'}; 
                    yl{1} = {'Total ocean';'reservoirs (mol N)'};                       yl{2} = {'Total ocean';'fluxes (mol N/yr)'}; 

                case 'Silcyc'
                    dim = [3 2];
                    lgpos = 'south'; ubound = [5e13 1e14 1e13 5e13 1e14 5e13];%[5e13 1e14 1e13 5e12 1e13 4e12 10]; 
                    lbound(5) = 1e12;                                                   lbound(2) = 1e11;%1e10;
                    xv{5} = {flux.wthr.sil,flux.sfw,tf.revweather};                     xv{2} = {flux.wthr.carb,tf.burial.CaCO3};  
                    ref{5}= {fluxx.wthr.sil,fluxx.t.sfw,fluxx.t.revweather};            ref{2}= {fluxx.wthr.carb, fluxx.burial.CaCO3.z+fluxx.burial.CaCO3.n}; 
                    rng{5}= {rflux.wthr.sil,rflux.t.sfw,rflux.t.revweather};            rng{2}= {rflux.wthr.carb, rflux.burial.CaCO3.z+rflux.burial.CaCO3.n}; 
                    c{5}  = {v.color.CO2,v.color.CaCO3,v.color.TA};                     c{2}  = {v.color.CaCO3, v.color.C};
                    nm{5} = {'subaerial','seafloor','reverse'};                         nm{2} = {'weathering','burial'}; 
                    yl{5} = {'Silicate weathering';'fluxes (mol C/yr)'};                yl{2} = {'CaCO_3';'fluxes (mol C/yr)'}; 

                    lbound(1) = 1e7;                                                    lbound(4) = 1e7; 
                    xv{1} = {flux.wthr.oxi,tf.burial.OC};                               xv{4} = {flux.volc.carb+flux.volc.CaCO3, flux.volc.OC};  
                    ref{1}= {fluxx.wthr.oxi,fluxx.burial.OC.z+fluxx.burial.OC.n};       ref{4}= {fluxx.volc.carb, fluxx.volc.OC}; 
                    rng{1}= {rflux.wthr.oxi,rflux.burial.OC.z+rflux.burial.OC.n};       rng{4}= {rflux.volc.carb, rflux.volc.OC}; 
                    c{1}  = {v.color.OC, v.color.C};                                    c{4}  = {v.color.CaCO3, v.color.OC};
                    nm{1} = {'weathering','burial'};                                    nm{4} = {'CaCO_3','Org C'}; 
                    yl{1} = {'Organic C';'fluxes (mol C/yr)'};                          yl{4} = {'Volcanism';'fluxes (mol C/yr)'}; 

                    lbound(3) = 1e8;                                                    lbound(6) = 1e8; 
                    xv{3} = {flux.meta.CaCO3,flux.meta.OC};                             xv{6} = {flux.mantle.CO2,flux.mantle.CH4};
                    ref{3}= {fluxx.meta.CaCO3,fluxx.meta.OC};                           ref{6}= {cfp.mantle.t.CO2.wa12, cfp.mantle.t.CH4.b04}; 
                    rng{3}= {rflux.meta.CaCO3,rflux.meta.OC};                           rng{6}= {cfr.mantle.t.CO2.wa12{1}, cfr.mantle.t.CH4.b04{1}};
                    c{3}  = {v.color.CaCO3, v.color.OC};                                c{6}  = {v.color.CO2,v.color.CH4};  
                    nm{3} = {'CaCO_3','Org C'};                                         nm{6} = {'CO_2','CH_4'}; 
                    yl{3} = {'Metamorphism';'fluxes (mol C/yr)'};                       yl{6} = {'Mantle';'outgassing (mol C/yr)'};

                case 'geoOC'
                    dim = [5 2]; ubound = [1e17 5e14 1e22 5e13 1e14 1e16 1e16 5e13 1e21 1e13]; lbound = [5e12 1e7 1e14 1e5 1e10 1e9 1e12 1e4 1e14 1e5]; lgpos = 'south';
                    xv{1} = {r.a.CH4};                                                  xv{2} = {gasex.CH4, flux.methox, flux.Hesc.CH4, flux.mantle.CH4};  
                    ref{1}= {v.atm.CH4pal};                                             ref{2}= {pref, fluxx.t.methox.*0.5, fluxx.t.Hesc.*4, cfp.mantle.t.CH4.b04}; 
                    rng{1}= {prng};                                                     rng{2}= {prng, rflux.t.methox.*0.5, rflux.t.Hesc.*4, cfr.mantle.t.CH4.b04{1}}; 
                    c{1}  = {v.color.CH4};                                              c{2}  = {v.color.CH4, v.color.O2, v.color.FeO, v.color.CH4}; 
                    lSty{1} = {'-'};                                                    lSty{2} = {'-','-','-',':'};
                    nm{1} = {''};                                                       nm{2} = {'air-sea','methox','H escape','mantle'};
                    yl{1} = {'Atmospheric';'reservoir (mol CH_4)'};                     yl{2} = {'Atmospheric';'fluxes (mol CH_4/yr)'}; 

                    xv{3} = {r.c.OC};                                                   xv{4} = {flux.wthr.oxi, flux.meta.OC};  
                    ref{3}= {rx.c.OC};                                                  ref{4}= {fluxx.wthr.oxi, fluxx.meta.OC}; 
                    rng{3}= {rxr.c.OC};                                                 rng{4}= {rflux.wthr.oxi, rflux.meta.OC}; 
                    c{3}  = {v.color.OC};                                               c{4}  = {v.color.OC,v.color.OC}; 
                    lSty{3} = {'-'};                                                    lSty{4} = {lsty.s,lsty.o};
                    nm{3} = {''};                                                       nm{4} = {'weathering','metamorphism'};
                    yl{3} = {'Continental organic';'carbon (mol C)'};                   yl{4} = {'Continental';'fluxes (mol C/yr)'}; 

                    xv{5} = {r.s.OC,r.d.OC};                                            xv{6} = {flux.export.OC, flux.sed.OC.n, flux.sed.OC.z};
                    ref{5}= {rx.s.OC,rx.d.OC};                                          ref{6}= {fluxx.export.OC,fluxx.sed.OC.n,fluxx.sed.OC.z}; 
                    rng{5}= {rxr.s.OC,rxr.d.OC};                                        rng{6}= {rflux.export.OC,rflux.sed.OC.n,rflux.sed.OC.z}; 
                    c{5}  = {v.color.OC,v.color.OC};                                    c{6}  = {v.color.OC,v.color.OC,v.color.OC};
                    lSty{5}= {lsty.s,lsty.d};                                           lSty{6} = {lsty.s,lsty.n,lsty.z};
                    nm{5} = {'surface','deep'};                                         nm{6} = {'export','nSed','zSed'}; 
                    yl{5} = {'Ocean organic';'carbon reservoirs (mol C)'};              yl{6} = {'Ocean';' fluxes (mol C/yr)'};

                    xv{7} = {r.n.OC,r.z.OC};                                            xv{8} = {flux.metha.OC.n+flux.denit.OC.n, flux.ammon.OC.n, flux.metha.OC.z+flux.denit.OC.z, flux.ammon.OC.z, flux.burial.OC.n, flux.burial.OC.z};
                    ref{7}= {pref,pref};                                                ref{8}= {pref,pref,pref,pref,fluxx.burial.OC.n,fluxx.burial.OC.z};
                    rng{7}= {rxr.n.OC,rxr.z.OC};                                        rng{8}= {prng,prng,prng,prng,rflux.burial.OC.n,rflux.burial.OC.z};
                    c{7}  = {v.color.OC,v.color.OC};                                    c{8}  = {v.color.C,v.color.CO2,v.color.C,v.color.CO2,v.color.OC,v.color.OC};
                    lSty{7}= {lsty.n,lsty.z};                                           lSty{8} = {lsty.n,lsty.n,lsty.z,lsty.z,lsty.n,lsty.z};
                    nm{7} = {'nOC','zOC'};                                              nm{8} = {'nAnox','nOxic','zAnox','zOxic','nBur','zBur'}; 
                    yl{7} = {'Reactive Sediment';'reservoirs (mol C)'};                 yl{8} = {'Reactive sediment';'fluxes (mol C/yr)'};

                    xv{9} = {r.u.OC};                                                   xv{10} = {flux.acc.OC, flux.subduct.OC, flux.volc.OC};
                    ref{9}= {rx.u.OC};                                                  ref{10}= {fluxx.acc.OC, fluxx.subduct.OC, fluxx.volc.OC};
                    rng{9}= {rxr.u.OC};                                                 rng{10}= {rflux.acc.OC, rflux.subduct.OC, rflux.volc.OC};
                    c{9}  = {v.color.OC};                                               c{10}  = {v.color.OC,v.color.OC,v.color.OC};
                    lSty{9} = {lsty.u};                                                 lSty{10} = {lsty.c,lsty.o,lsty.u}; 
                    nm{9} = nm{3};                                                      nm{10} = {'accretion','subduction','volcanism'};
                    yl{9} = {'Abyssal sediment';'organic carbon (mol C)'};              yl{10} = {'Subduction zone';'fluxes (mol C/yr)'}; 

                case 'geocarb'
                    dim = [5 2]; ubound = [1e19 5e14 5e22 1e14 1e20 5e14 1e17 1e14 5e21 2e13]; lbound = [1e16 1e10 1e20 1e10 1e10 5e11 1e14 1e10 1e18 5e10]; lgpos = 'south';
                    xv{1} = {r.a.CO2};                                                  xv{2} = {gasex.CO2, flux.methox, flux.mantle.CO2};  
                    ref{1}= {v.atm.CO2pal};                                             ref{2}= {pref, fluxx.t.methox.*0.5, cfp.mantle.t.CO2.wa12}; 
                    rng{1}= {prng};                                                     rng{2}= {prng, rflux.t.methox.*0.5, cfr.mantle.t.CO2.wa12{1}}; 
                    c{1}  = {v.color.CO2};                                              c{2}  = {v.color.CO2, v.color.O2, v.color.CO2}; 
                    lSty{1} = {'-'};                                                    lSty{2} = {'-','-',':'};
                    nm{1} = {''};                                                       nm{2} = {'air-sea','methox','mantle'};
                    yl{1} = {'Atmospheric';'reservoir (mol CO_2)'};                     yl{2} = {'Atmospheric';'fluxes (mol CO_2/yr)'}; 

                    xv{3} = {r.c.CaCO3, r.c.SiO3};                                      xv{4} = {flux.wthr.sil, flux.wthr.carb, flux.meta.CaCO3};  
                    ref{3}= {rx.c.CaCO3, v.ea.Si};                                      ref{4}= {fluxx.wthr.sil, fluxx.wthr.carb, fluxx.meta.CaCO3}; 
                    rng{3}= {rxr.c.CaCO3, prng};                                        rng{4}= {rflux.wthr.sil, rflux.wthr.carb, rflux.meta.CaCO3}; 
                    c{3}  = {v.color.CaCO3, v.color.CO2};                               c{4}  = {v.color.CO2, v.color.CaCO3,v.color.CaCO3}; 
                    lSty{3} = {'-','-'};                                                lSty{4} = {lsty.s,lsty.s,lsty.o};
                    nm{3} = {'CaCO_3','SiO_3'};                                         nm{4} = {'wSil','wCaCO_3','metamorphism'};
                    yl{3} = {'Continental';'reservoirs (mol)'};                         yl{4} = {'Continental';'fluxes (mol C/yr)'}; 

                    xv{5} = {r.s.CaCO3, r.s.DIC, r.d.CaCO3, r.d.DIC};                   xv{6} = {tf.precip, tf.diss, flux.sed.CaCO3.n, flux.sed.CaCO3.z};
                    ref{5}= {rx.s.CaCO3, rx.s.DIC, rx.d.CaCO3, rx.d.DIC};               ref{6}= {pref, pref, pref, pref}; 
                    rng{5}= {rxr.s.CaCO3, rxr.s.DIC, rxr.d.CaCO3, rxr.d.DIC};           rng{6}= {rflux.t.precip, rflux.t.diss, rflux.sed.CaCO3.n, rflux.sed.CaCO3.z}; 
                    c{5}  = {v.color.CaCO3,v.color.DIC, v.color.CaCO3,v.color.DIC};     c{6}  = {v.color.DIC,v.color.TA,v.color.CaCO3, v.color.CaCO3}; 
                    lSty{5}= {lsty.s,lsty.s,lsty.d,lsty.d};                             lSty{6} = {lsty.s,lsty.d,lsty.n,lsty.z};
                    nm{5} = {'sCaCO_3','sDIC','dCaCO_3','dDIC'};                        nm{6} = {'precip','diss','nSed','zSed'}; 
                    yl{5} = {'Ocean';'reservoirs (mol C)'};                             yl{6} = {'Ocean';' fluxes (mol C/yr)'};

                    xv{7} = {r.n.CaCO3,r.z.CaCO3};                                      xv{8} = {flux.burial.CaCO3.n,flux.diss.n, flux.burial.CaCO3.z,flux.diss.z,flux.sfw, tf.revweather};
                    ref{7}= {rx.n.CaCO3,rx.z.CaCO3};                                    ref{8}= {fluxx.burial.CaCO3.n,pref,fluxx.burial.CaCO3.z,pref, 5e12, fluxx.t.revweather};
                    rng{7}= {rxr.n.CaCO3,rxr.z.CaCO3};                                  rng{8}= {rflux.burial.CaCO3.n,prng,rflux.burial.CaCO3.z,prng, prng, rflux.t.revweather};
                    c{7}  = {v.color.CaCO3,v.color.CaCO3};                              c{8}  = {v.color.CaCO3,v.color.CO2,v.color.CaCO3,v.color.CO2, v.color.TA,v.color.TA};
                    lSty{7}= {lsty.n,lsty.z};                                           lSty{8} = {lsty.n,lsty.n,lsty.z,lsty.z,lsty.o,lsty.z};
                    nm{7} = {'nCaCO_3','zCaCO_3'};                                      nm{8} = {'nBur','nDiss','zBur','zDiss','SFW','RW'}; 
                    yl{7} = {'Reactive Sediment';'reservoirs (mol C)'};                 yl{8} = {'Reactive sediment';'fluxes (mol C/yr)'};

                    xv{9} = {r.u.CaCO3,r.o.CaCO3};                                      xv{10} = {flux.acc.CaCO3+flux.acc.carb, flux.subduct.CaCO3+flux.subduct.carb, flux.volc.CaCO3+flux.volc.carb};
                    ref{9}= {rx.u.CaCO3,rx.o.CaCO3};                                    ref{10}= {fluxx.acc.CaCO3, fluxx.subduct.CaCO3, fluxx.volc.carb};
                    rng{9}= {rxr.u.CaCO3,rxr.o.CaCO3};                                  rng{10}= {rflux.acc.CaCO3, rflux.subduct.CaCO3, rflux.volc.carb};
                    c{9}  = {v.color.CaCO3,v.color.CaCO3};                              c{10}  = {v.color.CaCO3,v.color.CaCO3,v.color.CO2};
                    lSty{9} = {lsty.u,lsty.o};                                          lSty{10} = {lsty.c,lsty.o,lsty.u}; 
                    nm{9} = {'uCaCO_3','oCaCO_3'};                                      nm{10} = {'accretion','subduction','volcanism'};
                    yl{9} = {'Abyssal sediments and';'slab reservoirs (mol C)'};        yl{10} = {'Subduction zone';'fluxes (mol C/yr)'}; 

                case 'geoN'
                    dim = [5 2]; ubound = [1e22 1e13 5e20 5e12 1e17 5e14 1e15 5e12 1e20 1e12]; lbound = [1e6 1e0 1e15 1e5 1e6 1e8 1e6 1e6 1e14 1e2]; lgpos = 'south'; 
                    xv{1} = {r.a.N2, r.a.NH3};                                          xv{2} = {gasex.N2, gasex.NH3, flux.pholys, flux.ammox, flux.Hesc.NH3, flux.mantle.N};  
                    ref{1}= {v.atm.N2pal,v.atm.NH3pal};                                 ref{2}= {pref, pref, pref, pref, pref, pref};
                    rng{1}= {prng,prng};                                                rng{2}= {prng, prng, prng, prng, prng, prng}; 
                    c{1}  = {v.color.N2,v.color.NH3};                                   c{2}  = {v.color.N2, v.color.NH3,v.color.NH3, v.color.O2, v.color.FeO, v.color.N}; 
                    lSty{1} = {'-','-'};                                                lSty{2} = {'-','-',':','-','-','--'};
                    nm{1} = {'N_2','NH_3'};                                             nm{2} = {'air-sea N_2','air-sea NH_3','pholys','ammox','H escape','mantle'};
                    yl{1} = {'Atmospheric';'reservoir (mol N)'};                    	yl{2} = {'Atmospheric';'fluxes (mol/yr)'}; 

                    xv{3} = {r.c.NH4, r.c.ON};                                          xv{4} = {flux.wthr.NH4,flux.wthr.ON, flux.meta.NH4, flux.meta.ON}; 
                    ref{3}= {rx.c.NH4,rx.c.ON};                                         ref{4}= {fluxx.wthr.NH4, fluxx.wthr.ON, fluxx.meta.NH4, pref}; 
                    rng{3}= {rxr.c.NH4,rxr.c.ON};                                       rng{4}= {rflux.wthr.NH4, rflux.wthr.ON, rflux.meta.NH4, prng}; 
                    c{3}  = {v.color.NH4,v.color.ON};                                   c{4}  = {v.color.NH4,v.color.ON,v.color.NH4,v.color.ON};
                    lSty{3} = {'-','-'};                                                lSty{4} = {lsty.s,lsty.s,lsty.o,lsty.o}; 
                    nm{3} = {'NH_4','ON'};                                              nm{4} = {'wNH_4','wON','mNH_4','mON'};
                    yl{3} = {'Continental';'reservoirs (mol N)'};                       yl{4} = {'Continental';'fluxes (mol N/yr)'}; 

                    xv{5} = {r.s.ON,r.s.RN, r.s.HNO3,r.d.ON,r.d.RN, r.d.HNO3};          xv{6} = {flux.export.ON, flux.sed.ON.n,flux.sed.ON.z};
                    ref{5}= {pref, rx.s.RN, rx.s.HNO3,pref, rx.d.RN, rx.d.HNO3};        ref{6}= {fluxx.export.OC./v.const.CNratio,pref,pref}; 
                    rng{5}= {prng,rxr.s.RN,rxr.s.HNO3,prng,rxr.d.RN,rxr.d.HNO3};        rng{6}= {rflux.export.OC./v.const.CNratio,prng,prng};
                    c{5}  = {v.color.ON,v.color.RN,v.color.HNO3,v.color.ON,v.color.RN,v.color.HNO3};c{6}  = {v.color.N,v.color.ON,v.color.ON};  
                    lSty{5}= {lsty.s,lsty.s,lsty.s,lsty.d,lsty.d,lsty.d};               lSty{6} = {lsty.s,lsty.n,lsty.z}; 
                    nm{5} = {'sON','sNH_3 + NH_4^+','sHNO_3','dON','dNH_3 + NH_4^+','dHNO_3'};nm{6} = {'export','nSed','zSed'}; 
                    yl{5} = {'Ocean';'reservoirs (mol N)'};                             yl{6} = {'Ocean';' fluxes (mol N/yr)'};

                    xv{7} = {r.n.ON, r.n.RN, r.z.ON, r.z.RN};                           xv{8} = {flux.burial.ON.n,flux.burial.ON.z, flux.hyd};
                    ref{7}= {pref,pref,pref,pref};                                      ref{8}= {fluxx.burial.ON.n, fluxx.burial.ON.z, pref};
                    rng{7}= {prng,prng,prng,prng};                                      rng{8}= {rflux.burial.ON.n, rflux.burial.ON.z,prng};
                    c{7}  = {v.color.ON,v.color.NH4,v.color.ON,v.color.NH4};            c{8}  = {v.color.ON,v.color.ON,v.color.NH4};
                    lSty{7}= {lsty.n,lsty.n,lsty.z,lsty.z};                             lSty{8} = {lsty.n,lsty.z,lsty.o};
                    nm{7} ={'nON','nRN','zON','zRN'};                                   nm{8} = {'nON','zON','hyd'};
                    yl{7} = {'Reactive Sediment';'reservoirs (mol N)'};                 yl{8} = {'Reactive sediment';'fluxes (mol N/yr)'};

                    xv{9} = {r.o.NH4, r.u.ON};                                          xv{10} = {flux.acc.ON, flux.acc.NH4, flux.subduct.ON, flux.subduct.NH4, flux.volc.ON, flux.volc.NH4, flux.cryst.NH4+flux.cryst.ON};
                    ref{9}= {rx.o.NH4,rx.u.ON};                                         ref{10}= {fluxx.acc.ON, fluxx.acc.NH4, fluxx.subduct.ON,fluxx.subduct.NH4, fluxx.volc.ON,fluxx.volc.NH4, fluxx.cryst.NH4}; 
                    rng{9}= {rxr.o.NH4,rxr.u.ON};                                       rng{10}= {rflux.acc.ON, rflux.acc.NH4, rflux.subduct.ON,rflux.subduct.NH4, rflux.volc.ON,rflux.volc.NH4, rflux.cryst.NH4};
                    c{9}  = {v.color.NH4,v.color.ON};                                   c{10}  = {v.color.ON, v.color.NH4, v.color.ON, v.color.NH4, v.color.ON, v.color.NH4, v.color.NH4};
                    lSty{9} = {lsty.o,lsty.u};                                          lSty{10} = {lsty.c,lsty.c,lsty.o,lsty.o,lsty.u,lsty.u,lsty.n}; 
                    nm{9} = {'oNH_4','uON'};                                            nm{10} = {'uON acc','oNH_4 acc','uON sub','oNH_4 sub','uON volc','oNH_4 volc','oNH_4 cryst',};
                    yl{9} = {'Abyssal sediment and';'slab reservoirs (mol N)'};         yl{10} = {'Subduction zone';'fluxes (mol N/yr)'}; 

                case 'geoP'
                    dim = [4 2];  ubound = [5e20 5e11 1e16 1e13 1e14 1e12 1e19 1e11]; lbound = [1e14 1e4 5e7 1e5 1e7 1e4 1e12 1e2]; lgpos = 'south';
                    xv{1} = {r.c.SP, r.c.CP, r.c.OP};                                   xv{2} = {flux.wthr.SP,flux.wthr.CP,flux.wthr.OP,flux.meta.OP};  
                    ref{1}= {rx.c.SP,rx.c.CP, rx.c.OP};                                 ref{2}= {fluxx.wthr.SP,fluxx.wthr.CP,fluxx.wthr.OP,pref}; 
                    rng{1}= {rxr.c.SP,rxr.c.CP, rxr.c.OP};                              rng{2}= {rflux.wthr.SP,rflux.wthr.CP,rflux.wthr.OP,prng}; 
                    c{1}  = {v.color.SP,v.color.CP,v.color.OP};                         c{2}  = {v.color.SP,v.color.CP,v.color.OP,v.color.OP}; 
                    lSty{1} = {'-','-','-'};                                            lSty{2} = {lsty.s,lsty.s,lsty.s,lsty.o};
                    nm{1} = {'silicate P','carbonate P','organic P'};                   nm{2} = {'wSP','wCP','wOP','mOP'};
                    yl{1} = {'Continental';'reservoirs (mol P)'};                       yl{2} = {'Continental';'fluxes (mol P/yr)'}; 

                    xv{3} = {r.s.OP, r.s.H3PO4,r.d.OP, r.d.H3PO4};                      xv{4} = {flux.export.OP, flux.sorb.n, flux.sed.OP.n,flux.sorb.z, flux.sed.OP.z};
                    ref{3}= {pref, rx.s.H3PO4,pref, rx.d.H3PO4};                        ref{4}= {fluxx.export.OC./v.const.CPratio,pref,pref,pref,pref}; 
                    rng{3}= {prng, rxr.s.H3PO4,prng, rxr.d.H3PO4};                      rng{4}= {rflux.export.OC./v.const.CPratio,prng,prng,prng,prng}; 
                    c{3}  = {v.color.OP,v.color.CP,v.color.OP,v.color.CP};              c{4}  = {v.color.P, v.color.SP,v.color.OP,v.color.SP,v.color.OP}; 
                    lSty{3}= {lsty.s,lsty.s,lsty.d,lsty.d};                             lSty{4} = {'-',lsty.n,lsty.n,lsty.z,lsty.z};
                    nm{3} = {'sOP','sH_3PO_4','dOP','dH_3PO_4'};                        nm{4} = {'export OP','nSorb','nSed','zSorb','zSed'}; 
                    yl{3} = {'Ocean';'reservoirs (mol P)'};                             yl{4} = {'Ocean';'fluxes (mol P/yr)'};

                    xv{5} = {r.n.FePO4,r.n.OP,r.z.FePO4,r.z.OP};                        xv{6} = {flux.burial.FePO4.n, flux.burial.PO4.n, flux.burial.OP.n, flux.burial.FePO4.z, flux.burial.PO4.z, flux.burial.OP.z};
                    ref{5}= {pref,pref,pref,pref};                                      ref{6}= {pref,fluxx.burial.PO4.n,fluxx.burial.OP.n,pref,pref,fluxx.burial.OP.z};
                    rng{5}= {prng,prng,prng,prng};                                      rng{6}= {prng,rflux.burial.PO4.n,rflux.burial.OP.n,rflux.t.sorb,rflux.burial.CP.z,rflux.burial.OP.z};
                    c{5}  = {v.color.SP,v.color.OP,v.color.SP,v.color.OP};              c{6}  = {v.color.SP,v.color.CP,v.color.OP,v.color.SP,v.color.CP,v.color.OP}; 
                    lSty{5}= {lsty.n,lsty.n,lsty.z,lsty.z};                             lSty{6} = {lsty.n,lsty.n,lsty.n,lsty.z,lsty.z,lsty.z};
                    nm{5} = {'nFe-P','nOP','zFe-P','zOP'};                              nm{6} = {'nSP','nCP','nOP','zSP','zCP','zOP'}; 
                    yl{5} = {'Reactive Sediment';'reservoirs (mol P)'};                 yl{6} = {'Reactive sediment';'burial fluxes (mol P/yr)'};

                    xv{7} = {r.u.SP,r.u.CP,r.u.OP};                                     xv{8} = {flux.acc.OP, flux.acc.SP,flux.acc.CP, flux.subduct.OP, flux.subduct.SP,flux.subduct.CP, flux.cryst.OP+flux.cryst.SP+flux.cryst.CP};
                    ref{7}= {pref,pref, rx.u.OP};                                       ref{8}= {fluxx.acc.OP, fluxx.acc.SP,fluxx.acc.CP, fluxx.subduct.OP,fluxx.subduct.SP,fluxx.subduct.CP, fluxx.cryst.CP}; 
                    rng{7}= {prng,prng, rxr.u.OP};                                      rng{8}= {rflux.acc.OP, rflux.acc.SP,rflux.acc.CP,rflux.subduct.OP,rflux.subduct.SP,rflux.subduct.CP,rflux.cryst.CP};
                    c{7}  = {v.color.SP,v.color.CP,v.color.OP};                         c{8}  = {v.color.OP, v.color.SP, v.color.CP, v.color.OP, v.color.SP, v.color.CP, v.color.SP};
                    lSty{7} = {lsty.u,lsty.u,lsty.u};                                   lSty{8} = {lsty.u,lsty.u,lsty.u,lsty.o,lsty.o,lsty.o,lsty.c}; 
                    nm{7} = nm{1};                                                      nm{8} = {'OP acc','SP acc','CP acc','OP sub','SP sub','CP sub','SP cryst'};
                    yl{7} = {'Abyssal sediment';'reservoirs (mol P)'};                  yl{8} = {'Subduction zone';'fluxes (mol P/yr)'}; 

                case 'mods'
                    lbound(1) = 1e2; lbound(2) = 0; dim = 'v';  
                    xv{1} = r.a.CO2 .* 1e6 ./ v.atm.mol;                                xv{2} = {msil, mcarb, msfw./(msfw + Xsp)}; 
                    ref{1}= 400;                                                        ref{2}= {1, 1, 1}; 
                    rng{1}= [280 500];                                                  rng{2}= {prng,prng,prng};
                    c{1}  = '';                                                         c{2}  = {v.color.CO2, v.color.CaCO3, v.color.DIC};
                    nm{2} = '';                                                         nm{2} = {'silicate','carbonate','seafloor'};
                    yl{2} = {'Atmospheric ';'CO_2 (ppm)'};                              yl{2} = {'Weathering';'modifiers (S_x)'}; 

                case 'forcings'
                    lgpos = 'north'; dim = 'v'; frames = [1e8 4.25e9]; % shrink reference zone
                    lbound(2) = 0; ubound(2) = 1.2;                                     lbound(1) = 1e11; ubound(1) = 5e12;   
                    xv{2} = (v.ea.Si./r.c.SiO3).^-1.5;                                  xv{1} = flux.mantle.FeO; 
                    ref{2}= 1;                                                          ref{1}= 3e11; 
                    rng{2}= prng;                                                       rng{1}= [1 5].*1e11;
                    lSty{2} = {'-'};                                                    lSty{1} = {'-'}; 
                    c{2}  = '';                                                         c{1}  = ''; 
                    nm{2} = '';                                                         nm{1} = '';
                    yl{2} = {'Relative continental';'erosion (\chi_{eros})'};           yl{1} = {'Reductant';'outgassing (F_{mantle,FeO})'}; 

                    lbound(4) = 0; ubound(4) = 1;                                       lbound(3) = 1e3;  ubound(3) = 1400;   
                    xv{4} = {tdep.photo, tdep.fungi, tdep.sink, tdep.plant, tdep.terr}; xv{3} = rf.sc;  
                    ref{4}= {pref, pref, pref, pref, pref};                             ref{3}= v.S.Pref; 
                    rng{4}= {prng, prng, prng, prng, prng};                             rng{3}= [1360 1362]; 
                    lSty{4} = {'-','-','-','-',':'};                                    lSty{3} = {'-'}; 
                    c{4}  = {v.color.O2,v.color.SP,v.color.OC,v.color.CO2,v.color.T};   c{3}  = v.color.Modern; 
                    nm{4} = {'Photosynthesizers','Fungi','Body size','Land plants','Terrestrial productivity'};nm{3} = '';
                    yl{4} = {'Time Dependent';'forcings (\chi_{x}(t))'};                yl{3} = {'Solar Constant';'(W/m^2)'}; 


                case 'osys' % oxygen balance, gasex of oxygen and reductants, h escape, weathering/oxidation fluxes, etc. 
                    lgpos = 'southeast'; dim = [5 2]; nmrg = [9 10]; 
                    lbound = [1e-8 1e-10 1e2 1e8 1e9 1e10 1e6 1e9 1e10]; ubound = [1 1e5 1e15 1e12 5e16 5e16 1e13 1e13 5e16]; 
                    N = flux.prod.O2 - tf.ammon; red = flux.mantle.FeO./4; dC = tf.burial.OC - flux.wthr.oxi; % OC 1:1 and FeO 4:1 w/r/to oxygen draw
                    sources = flux.prod.O2 + flux.Hesc.H./4; sinks = tf.ammon + 2.*(tf.nitr + tf.mtrophy + flux.methox) + flux.wthr.oxi + 0.5.*flux.wthr.Feox + 0.75.*flux.ammox; 
                    xv{1} = {r.a.O2./v.atm.O2pal};                                     xv{2} = {(N ./ (red - dC)).^0.7}; 
                    ref{1}= {pref};                                                    ref{2}= {pref}; 
                    rng{1}= {prng};                                                    rng{2}= {prng};
                    c{1}  = {v.color.O2};                                              c{2}  = {''}; 
                    nm{1} = {''};                                                      nm{2} = {''};
                    lSty{1} = {'-'};                                                   lSty{2} = lSty{1};
                    yl{1} = {'Present Atmospheric';'level (mol O_2/mol O_2_{ref})'};   yl{2} = {'(N/(r - \DeltaC))^{0.7}'};

                    xv{3} = {flux.ammox.*0.75, flux.methox.*2};                         xv{4} = {flux.Hesc.H./4}; 
                    ref{3}= {pref,fluxx.t.methox.*2};                                   ref{4}= {fluxx.t.Hesc./4}; 
                    rng{3}= {prng,rflux.t.methox.*2};                                   rng{4}= {rflux.t.Hesc./4}; 
                    c{3}  = {v.color.NH3,v.color.CH4};                                  c{4}  = c{1};
                    nm{3} = {'NH_3','CH_4'};                                            nm{4} = nm{1};
                    lSty{3} = {'-','-'};                                                lSty{4} = lSty{1};
                    yl{3} = {'Photochemical oxidation';'fluxes (mol O_2 eq/yr)'};       yl{4} = {'Hydrogen Escape';'flux (mol O_2 eq/yr)'}; 

                    xv{5} = {tf.ammon, tf.nitr.*2, tf.mtrophy.*2,flux.photox.O2};       xv{6} = {flux.prod.O2}; 
                    ref{5}= {fluxx.t.oxrm,fluxx.t.nitr.*2,pref,pref};                   ref{6}= {fluxx.prod.O2}; 
                    rng{5}= {rflux.t.oxrm,rflux.t.nitr.*2,prng,prng};                   rng{6}= {rflux.prod.O2}; 
                    c{5}  = {v.color.OC,v.color.HNO3,v.color.CH4,v.color.Fe};           c{6}  = c{1};
                    nm{5} = {'ammonification','nitrification','methanotrophy','Fe photo-oxidation'};nm{6} = nm{1};
                    lSty{5} = {'-','-','-','-'};                                        lSty{6} = lSty{1};
                    yl{5} = {'Total ocean oxygen';'sinks (mol O_2 eq/yr)'};             yl{6} = {'Photosynthesis';'flux (mol O_2 eq/yr)'}; 

                    xv{7} = {flux.wthr.oxi,flux.wthr.Feox.*0.5};                        xv{8} = {tf.burial.OC}; 
                    ref{7}= {fluxx.wthr.oxi,pref};                                      ref{8}= {fluxx.burial.OC.z}; 
                    rng{7}= {rflux.wthr.oxi,prng};                                      rng{8}= {rflux.burial.OC.z}; 
                    c{7}  = {v.color.OC,v.color.Fe};                                    c{8}  = {v.color.OC}; 
                    nm{7} = {'Org Carbon','Reduced Iron'};                              nm{8} = {''};
                    lSty{7} = {'-','-'};                                                lSty{8} = lSty{1};  
                    yl{7} = {'Oxidative weathering';'fluxes (mol O_2 eq/yr)'};          yl{8} = {'Organic carbon';'burial (mol O_2 eq/yr)'}; 

                    xv{9} = {sources,sinks};                       
                    ref{9}= {pref, pref};                 
                    rng{9}= {prng,prng};               
                    c{9}  = {'',''};          
                    nm{9} = {'sources','sinks'};   
                    lSty{9} = {'-',':'};   
                    yl{9} = {'Oxygen flux';'balance (mol O_2 eq/yr)'};  

                case 'balance'  
                    dim = [4 3]; lgpos = 'south'; % atm/cont, surface, deep fluxes
                    pagesz = [20 12]; % bigger page for debugging
                    lbound = [1e6 1e6 1e4 1e11 1e9 1e8 1e8 1e6 1e6 1e11 1e12 1e12]; ubound = [1e15 1e17 5e16 1e14 5e16 5e16 1e12 5e14 1e13 1e15 1e15 1e15];  
                    % OXYGEN
                    xv{1} = {gasex.O2,flux.methox.*2,flux.ammox.*0.75,flux.wthr.oxi};   xv{2} = {mix.O2.s,mix.diff.O2.n,flux.prod.O2,flux.ammon.OC.s,flux.nitr.s.*2,flux.mtrophy.s.*2,flux.photox.O2}; 
                    ref{1}= {pref,fluxx.t.methox.*2,pref,fluxx.wthr.oxi};               ref{2}= {pref,pref,fluxx.prod.O2,fluxx.ammon.OC.s,fluxx.t.nitr,pref,pref}; 
                    rng{1}= {prng,rflux.t.methox.*2,pref,rflux.wthr.oxi};               rng{2}= {prng,prng,rflux.prod.O2,rflux.ammon.OC.s,rflux.t.nitr,prng,prng}; 
                    c{1}  = {v.color.O2,v.color.CH4,v.color.NH3,v.color.C};             c{2}  = {v.color.O2,v.color.O2,v.color.CO2,v.color.C,v.color.HNO3,v.color.CH4,v.color.FeO}; 
                    lSty{1} = {'-',':',':',':'};                                        lSty{2} = {'-',':','-',':',':',':',':'};
                    nm{1} =  {'O_2 gas exchange','CH_4 oxidation','NH_3 oxidation','oxidative weathering'};nm{2} = {'mixing','diffusion','photosynthesis','ammonification','nitrification','methanotrophy','photo-oxidation'}; 
                    yl{1} = {'Atmospheric fluxes (mol O_2/yr)'};                        yl{2} = {'Surface ocean fluxes (mol O_2/yr)'}; 

                    xv{3} = {mix.O2.d,mix.diff.O2.z,flux.ammon.OC.d,flux.nitr.d.*2,flux.mtrophy.d.*2};
                    ref{3}= {pref,pref,fluxx.ammon.OC.d,pref,pref}; 
                    rng{3}= {prng,pref,rflux.ammon.OC.d,prng,prng}; 
                    c{3}  = {v.color.O2,v.color.O2,v.color.C,v.color.HNO3,v.color.CH4}; 
                    lSty{3} = {'-',':',':',':',':'};
                    nm{3} = {'mixing','diffusion','ammonification','nitrification','methanotrophy'}; 
                    yl{3} = {'Deep ocean fluxes (mol O_2/yr)'};
                    % ORGANIC CARBON
                    xv{4} = {gasex.CH4};                                                xv{5} = {flux.prod.CO2,flux.export.OC,flux.sed.OC.n,flux.ammon.OC.s+flux.denit.OC.s+flux.metha.OC.s};
                    ref{4}= {pref};                                                     ref{5}= {fluxx.prod.CO2,fluxx.export.OC,fluxx.sed.OC.n,fluxx.ammon.OC.s+(fluxx.t.denit.*(106./84.8))+fluxx.metha.OC.s};
                    rng{4}= {prng};                                                     rng{5}= {rflux.prod.CO2,rflux.export.OC,rflux.sed.OC.n,rflux.ammon.OC.s+(rflux.t.denit.*(106./84.8))+rflux.metha.OC.s};
                    c{4}  = {v.color.CH4};                                              c{5}  = {v.color.OC,v.color.OC,v.color.DIC,v.color.C}; 
                    lSty{4} = {'-'};                                                    lSty{5} = {'-',':',':',':'};
                    nm{4} = {''};                                                       nm{5} = {'photosynthesis','export','sed','remineralization'}; 
                    yl{4} = {'Air-sea exchange flux (mol CH_4/yr)'};                    yl{5} = {'Surface ocean fluxes (mol C/yr)'}; 

                    xv{6} = {flux.sed.OC.z,flux.ammon.OC.d+flux.denit.OC.d+flux.metha.OC.d};
                    ref{6}= {fluxx.sed.OC.z,fluxx.ammon.OC.d+(fluxx.t.denit.*(106./84.8))+fluxx.metha.OC.d}; 
                    rng{6}= {rflux.sed.OC.z,rflux.ammon.OC.d+(rflux.t.denit.*(106./84.8))+rflux.metha.OC.d};
                    c{6}  = {v.color.O2,v.color.OC}; 
                    lSty{6} = {':',':'};
                    nm{6} = {'sed','remineralization'}; 
                    yl{6} = {'Deep ocean fluxes (mol C/yr)'};
                    % PHOSPHORUS
                    xv{7} = {flux.wthr.SP,flux.wthr.CP,flux.wthr.OP};                   xv{8} = {mix.H3PO4.s,mix.diff.H3PO4.n,flux.prod.H3PO4,flux.export.OP,flux.sed.OP.n,(flux.ammon.H3PO4.s+flux.denit.H3PO4.s+flux.metha.H3PO4.s),flux.sorb.s+flux.sorb.n};
                    ref{7}= {fluxx.wthr.SP,fluxx.wthr.CP,fluxx.wthr.OP};                ref{8}= {pref,pref,fluxx.prod.CO2./v.const.CPratio,fluxx.export.OC./v.const.CPratio,fluxx.sed.OP.n,(fluxx.ammon.H3PO4.s+(fluxx.t.denit.*(16./84.8))+fluxx.metha.H3PO4.s),pref};
                    rng{7}= {rflux.wthr.SP,rflux.wthr.CP,rflux.wthr.OP};                rng{8}= {prng,prng,rflux.prod.CO2./v.const.CPratio,rflux.export.OC./v.const.CPratio,rflux.sed.OP.n,(rflux.ammon.H3PO4.s+(rflux.t.denit.*(16./84.8))+rflux.metha.H3PO4.s),prng}; 
                    c{7}  = {v.color.SP,v.color.CP,v.color.OP};                         c{8}  = {v.color.H3PO4,v.color.H3PO4,v.color.CO2,v.color.P,v.color.OP,v.color.C,v.color.SP}; 
                    lSty{7} = {'-','-','-'};                                            lSty{8} = {'-',':',':',':',':','-',':'};
                    nm{7} =  {'silicate P','carbonate P','organic P'};                  nm{8} = {'mixing','diffusion','photosynthesis','export','sed','remineralization','Fe-sorption'}; 
                    yl{7} = {'Continental weathering fluxes (mol P/yr)'};               yl{8} = {'Surface ocean fluxes (mol P/yr)'}; 

                    xv{9} = {mix.H3PO4.d,mix.diff.H3PO4.z,flux.sed.OP.z,(flux.ammon.H3PO4.d+flux.denit.H3PO4.d+flux.metha.H3PO4.d),flux.sorb.d+flux.sorb.z};
                    ref{9}= {pref,pref,fluxx.sed.OP.z,(fluxx.ammon.H3PO4.d+(fluxx.t.denit.*(16./84.8))+fluxx.metha.H3PO4.d),pref};
                    rng{9}= {prng,prng,rflux.sed.OP.z,(rflux.ammon.H3PO4.d+(rflux.t.denit.*(16./84.8))+rflux.metha.H3PO4.d),prng};
                    c{9}  = {v.color.H3PO4,v.color.H3PO4,v.color.OP,v.color.OC,v.color.SP};
                    lSty{9} = {'-',':',':',':',':'};
                    nm{9} = {'mixing','diffusion','sed','remineralization','Fe-sorption'}; 
                    yl{9} = {'Deep ocean fluxes (mol P/yr)'};
                    % CARBONATE
                    xv{10} = {gasex.CO2};                                               xv{11} = {flux.precip.s,flux.export.CaCO3,flux.sed.CaCO3.n};
                    ref{10}= {pref};                                                    ref{11}= {fluxx.t.precip,fluxx.export.CaCO3,fluxx.sed.CaCO3.n};
                    rng{10}= {prng};                                                    rng{11}= {rflux.t.precip,rflux.export.CaCO3,rflux.sed.CaCO3.n}; 
                    c{10}  = {v.color.CO2};                                             c{11}  = {v.color.CO3,v.color.CaCO3,v.color.C}; 
                    lSty{10} = {'-'};                                                   lSty{11} = {'-',':',':'};
                    nm{10} =  {'CO_2 gas exchange'};                                    nm{11} = {'precipitation','export','sed'}; 
                    yl{10} = {'Air-sea exchange flux (mol CO_2/yr)'};                   yl{11} = {'Surface ocean fluxes (mol P/yr)'}; 

                    xv{12} = {flux.diss.d+flux.diss.z,flux.sed.CaCO3.z};
                    ref{12}= {fluxx.t.diss,fluxx.sed.CaCO3.z};
                    rng{12}= {rflux.t.diss,rflux.sed.CaCO3.z};
                    c{12}  = {v.color.TA,v.color.CaCO3};
                    lSty{12} = {':',':'};
                    nm{12} = {'dissolution','sed'}; 
                    yl{12} = {'Deep ocean fluxes (mol P/yr)'};

                case 'mantlesurf' % carbon/nitrogen surface flux balance  and mantle reservoir changes
                    lgpos = 'south'; dim = 'v'; lbound = [1e-3 1e-5]; ubound = [1 1];
                    cC = r.c.CaCO3+r.c.OC+r.u.CaCO3+r.u.OC+r.o.CaCO3; oaC = tr.DIC+tr.CaCO3+tr.CH4+tr.OC+r.a.CO2+r.a.CH4+r.s.LB;
                    cN = r.c.NH4+r.c.ON+r.u.ON+r.o.NH4; oaN = tr.fixedN+tr.ON+2.*(tr.N2+r.a.N2)+r.a.NH3+(r.s.LB./v.const.CNratio);
                    rxTC =  rx.m.C+rx.c.C+rx.ao.C; rxTN = rx.m.N+rx.c.N+rx.ao.N; rxrTC = rxr.m.C+rxr.c.C+rxr.ao.C; rxrTN = rxr.m.N+rxr.c.N+rxr.ao.N; % totals for reservoir estimates
                    xv{2} = {r.m.C./totalres.C,cC./totalres.C,oaC./totalres.C};        xv{1} = {r.m.N./totalres.N,cN./totalres.N,oaN./totalres.N}; 
                    ref{2}= {rx.m.C./rxTC,rx.c.C./rxTC,rx.ao.C./rxTC};                 ref{1}= {rx.m.N./rxTN,rx.c.N./rxTN,rx.ao.N./rxTN}; 
                    rng{2}= {rxr.m.C./rxrTC,rxr.c.C./rxrTC,rxr.ao.C./rxrTC};           rng{1}= {rxr.m.N./rxrTN,rxr.c.N./rxrTN,rxr.ao.N./rxrTN};
                    c{2}  = {v.color.C,v.color.C,v.color.C};                           c{1}  = {v.color.N,v.color.N,v.color.N}; 
                    lSty{2} = {'-','--',':'};                                          lSty{1} = lSty{2};
                    nm{2} = {'mantle','crust','atmosphere-ocean'};                     nm{1} = nm{2};
                    yl{2} = {'Apportionment of total Carbon (mol C/total C)'};         yl{1} = {'Apportionment of total Nitrogen (mol N/total N)'};
            end
            len = length(xv);
            switch figfocus{ig}
                case {'temprf','osys','bio'}
                   len = len + 1;
                otherwise 
                    nmrg = 'x'; 
            end
            % generate subplots and change spacing
            [fg,ax] = plot_GenerateSubplots(len,dim,nmrg,'x'); 
            if ~ischar(dim) 
                plot_SqueezeSubplots(ax,'v',dim(2),'tight','n','y');
                refspace = 7e7; % more space between references, so as to not overlap
            elseif ischar(dim) 
                plot_SqueezeSubplots(ax,dim,1,'tight','n','y');
                refspace = 4e7; % less spacing required because of wider reference zone
            end
            for iv = 1:length(xv)
                if ~iscell(xv{iv}) % single line
                    if ischar(c{iv})                                    % undesignated colors
                       crd = get(gca,'ColorOrder');                     % get the color just used to apply to the reference
                       c{iv} = crd(iv,:); 
                    end
                    plot(ax(iv),t,xv{iv},'color',c{iv});
                    PlotRefRange(ax(iv),refx,ref{iv},rng{iv},c{iv});   % plot reference and range
                else % multiples lines
                   count = 1;
                   for ixv = 1:length(xv{iv})
                       if exist('lSty','var')                           % if assigned line styles, use them
                           LINE = lSty{iv}{ixv};
                       elseif strcmp(nm{iv}{ixv},'fixed N')
                           LINE = ':'; 
                       else
                           LINE = '-';
                       end
                       legname = nm{iv}{ixv};
                       if ischar(c{iv}{ixv})                            % for undesignated colors
                            crd = get(gca,'ColorOrder');                % get the color just used to apply to the reference
                            c{iv}{ixv} = crd(ixv,:); 
                       end
                       pn = plot(ax(iv),t,xv{iv}{ixv},LINE,'color',c{iv}{ixv}); 
                       if ~strcmp(nm{iv}{ixv},'')
                           pn.DisplayName = legname;
                       else 
                           pn.HandleVisibility = 'off'; 
                       end
                       if ixv == 1 && (ref{iv}{ixv} == pref) && (rng{iv}{ixv}(1) == prng(1)) % ensures that the real references are places closest to the line instead of pushed to the right
                           reftm = refx; 
                       elseif ixv > 1 && (ref{iv}{ixv} == pref) && (rng{iv}{ixv}(1) == prng(1)) % ensures that the real references are places closest to the line instead of pushed to the right
                           reftm = refx; 
                       else 
                           reftm = refx + count.*refspace;               % so the references plot side by side, not overlapping
                           count = count + 1;
                       end
                       PlotRefRange(ax(iv),reftm,ref{iv}{ixv},rng{iv}{ixv},c{iv}{ixv},LINE); % plot reference and range
                   end
                   children = get(ax(iv),'children');                   % if the subplot contains more than one line, produce a legend
                   if length(children) > 1  
                       if strcmp(children(1).DisplayName,'z')           % this is a box, only show this legend on the first subplot
                           if iv == 1
                               lg = legend(ax(iv),'-DynamicLegend','location',lgpos);
                               lg.NumColumns = 2;
                           end
                       else                                             % print a lengend on every subplot that has multiple lines
                           lg = legend(ax(iv),'-DynamicLegend','location',lgpos,'FontSize',fontsz-1);
                           if length(children) > 2 && length(children) <= 6 
                               lg.NumColumns = round(ixv/2); 
                           elseif length(children) > 6
                               lg.NumColumns = round(ixv/3); 
                           else
                               lg.NumColumns = length(children);  
                           end
                           lsizes = get(lg,'ItemTokenSize');            % get the sizes of line icons in the legend
                           lg.ItemTokenSize = lsizes./2;                % divide the size vector in half
                       end
                   end
                end
                % adjust the x and y scales according to dataset, and widen if needed
                oyl = get(ax(iv),'ylim'); 
                if oyl(2) - oyl(1) < 500                                % keep a lin scale for this range
                    if strcmp(figfocus{ig},'ocean') && iv ~= 10
                        set(ax(iv),'yscale','log'); oyl = get(ax(iv),'ylim'); 
                        SetYLim(ax(iv),oyl,lbound(iv),ubound(iv));
                    elseif strcmp(figfocus{ig},'temprf') && iv == 1
                        set(ax(iv),'yscale','log'); oyl = get(ax(iv),'ylim'); 
                        SetYLim(ax(iv),oyl,lbound(iv),ubound(iv));
                    elseif strcmp(figfocus{ig},'forcings') && iv == 2
                        set(ax(iv),'yscale','log'); oyl = get(ax(iv),'ylim'); 
                        SetYLim(ax(iv),oyl,lbound(iv),ubound(iv));
                    elseif strcmp(figfocus{ig},'osys') && iv == 1
                        set(ax(iv),'yscale','log'); oyl = get(ax(iv),'ylim'); 
                        SetYLim(ax(iv),oyl,lbound(iv),ubound(iv));
                    elseif strcmp(figfocus{ig},'bio') && iv >= 5
                        set(ax(iv),'yscale','log'); 
                    elseif strcmp(figfocus{ig},'mantlesurf') 
                        set(ax(iv),'yscale','log');
                    else
                        set(ax(iv),'yscale','lin');
                    end
                elseif strcmp(figfocus{ig},'ocean') && iv == 8
                    set(ax(iv),'yscale','lin');
                elseif strcmp(figfocus{ig},'forcings') && iv == 1
                    set(ax(iv),'yscale','lin');
                else  % set log scale
                    set(ax(iv),'yscale','log'); oyl = get(ax(iv),'ylim'); 
                    SetYLim(ax(iv),oyl,lbound(iv),ubound(iv)); 
                end
                set(ax(iv),'xscale','lin','xlim',frames,'FontSize',fontsz);
                ylabel(ax(iv),yl{iv},'FontSize',fontsz);                     
               if iv == 1                                               % only label the first subplot 
                  lblon = 'on'; rflbl = 'References';
               else
                  lblon = 'off'; rflbl = ''; 
               end
               if ~strcmp(figfocus{ig},'forcings')
                   plot_Transitions(ax(iv),v,xlinet,lblon,fontsz);
               else                                                     % show reference zone only
                   xline(ax(iv),xlinet,'-',rflbl,'LabelVerticalAlignment','middle',...
                       'LabelHorizontalAlignment','right','LabelOrientation','aligned',...
                       'FontSize',fontsz,'HandleVisibility','off'); 
               end
               if ~ischar(dim)                                          % assign xtick labels to only the final subplots, depending on number of columns
                   if iv >= length(xv)-(dim(2)-1) 
                        plot_Agescale(ax(iv),'ga');  
                        ax(iv).XLabel.FontSize = 10; 
                        ax(iv).XTickLabels{end} = '';
                   end
               else
                   if iv == length(xv)
                        plot_Agescale(ax(iv),'ga');
                        ax(iv).XLabel.FontSize = 10; 
                        ax(iv).XTickLabels{end} = '';
                   end
               end
            end
        end
        if length(ax) > 1                                               % add labels to each subplot
            if ~ischar(dim)
                plot_LabelSubplots(ax,'out','alpha',fontsz);
            else
                plot_LabelSubplots(ax,'out','alpha',fontsz,1);
            end
        end
        % Generate a PDF! NOTE: the folder you save to must already exist!
        set(fg,'Units','inches','PaperUnits','inches','PaperSize',pagesz);
        PrintPDFToFolder(pagesz(1),pagesz(2),['NominalRun',figfocus{ig}],v.figfolder);
    end
end

%% Subfunction: plot references 

% Inputs are axis handle (ax), time OR x value to plot at (xval), the
% reference value (ref), the range for the error bar (rng), symbol color
% (color) and a name to display (disp) or linestyle for error bar, if 
% desired (optional inputs). 
% Range input should obviously be a vector [low high]. 

function PlotRefRange(ax,xval,ref,rng,color,lsty,disp)
    if ~exist('lsty','var')
        lsty = '-'; 
    end
    if ~exist('disp','var')
       disp = 'n';  
    end
    yval = ref; % plot the reference
    if rng(1) == 1e-15 && ref ~= 1e-15% this is my "don't have a reference" placeholders
        yval = ref; 
        err = 'x';% don't plot error bar  
    elseif ref == 1e-15 && rng(1) ~= 1e-15 % if I only have a range, use the mean of that as single reference
        yval = mean(rng); 
        err = [yval-rng(1) rng(2)-yval]; % difference between the reference and the range values
    elseif ref ~= 1e-15 && rng(1) ~= 1e-15
        yval = ref; 
        err = [ref-rng(1) rng(2)-ref]; % difference between the reference and the range values
    else % plot nothing
        return 
    end
% plot the reference value as symbol corresponding with a line style and reference range as error bars
    switch lsty % different naming conventions THANKS MATLAB
        case '-'
            LSTYLE = 'solid'; mrk = 'o';
        case ':'
            LSTYLE = 'dotted';mrk = 's';
        case '--'
            LSTYLE = 'dashed';mrk = 'p';
        case '-.'
            LSTYLE = 'dashdot';mrk = '^';
    end   
    if ischar(err)
        eb = scatter(ax,xval,yval,30,color,mrk);
    else
        eb = errorbar(ax,xval,yval,err(1),err(2),mrk,'MarkerSize',5,'color',color);          
        eb.Bar.LineStyle = LSTYLE; 
    end
    switch disp
        case 'n'
            set(eb,'HandleVisibility','off');
        otherwise
            set(eb,'DisplayName',disp); 
    end
eb.LineWidth = 1.2;
end

%% Subfunction: print more descriptive species and flux names for legends/axis labels
% depending on the simple structure naming system, assign a descriptive
% string name to the flux for plotting in legends and axes labels

% SAMPLE USAGE:
%   enter a flux name and species :
%       MakeDescriptiveHandle('fix','N2')                    --> 'N_2 fixing'
%   enter a flux name, a reservoir, and a species
%       MakeDescriptiveHandle('denit',{'s','HNO3'})          --> 's HNO_3 denitrification remin'
%   enter a flux name, a reservoir, species, and an order
%       MakeDescriptiveHandle('denit',{'s','HNO3'},'rlast')  --> 'HNO_3 denitrification remin s'
%   enter a species and a reservoir
%        MakeDescriptiveHandle('SPECIES',{'s','HNO3'})       --> 's HNO_3'
%   enter ONLY a species
%        MakeDescriptiveHandle('SPECIES','HNO3')             --> ' HNO_3'
%   enter ONLY a flux
%        MakeDescriptiveHandle('nitr')                       --> ' NH_4 nitrification'

function fname = plot_DescriptiveHandle(FLUX,desig,order)
if ~exist('desig','var') || isempty(desig) % this is an optional input, but one may still want to include an order in a loop
    desig = 666; 
end
if ~exist('order','var') % also optional!
    order = 'rfirst';
end
% the first input is a flux name, the second input is either a species or reservoir (or both!)
switch FLUX
    case 'acc'
        FLUXNAME = 'accretion';
    case 'ammon'
        FLUXNAME = 'ammonification';% remin'; 
    case 'ammox'
        FLUXNAME = 'NH_3 photo-oxidation';
    case 'assim'
        FLUXNAME = 'photosynthesis';
    case 'boundary'
        FLUXNAME = 'boundary flux';
    case 'cryst'
        FLUXNAME = 'recrystallization';
    case 'denit'
        FLUXNAME = 'denitrification';% remin'; 
    case 'diff'
        FLUXNAME = 'diffusion'; 
    case 'diss'
        FLUXNAME = 'CaCO_3 dissolution';
    case 'ferrotrophy'
        FLUXNAME = 'photoferrotrophy';
    case {'fix','fixing'}
        FLUXNAME = 'N_2 fixing';  
    case 'forg'
        FLUXNAME = 'Org C burial fraction';
    case 'haze'
        FLUXNAME = 'aerosol haze production';
    case 'Hesc'
        FLUXNAME = 'hydrogen escape';
    case 'mantle'
        FLUXNAME = 'mantle reductant outgassing';
    case 'meta'
        FLUXNAME = 'metamorphic outgassing';
    case 'metha'
        FLUXNAME = 'methanogenesis';% remin';
    case 'methox'
        FLUXNAME = 'methane oxidation';
    case 'mtrophy'
        FLUXNAME = 'methanotrophy';
    case 'ncp'
        FLUXNAME = 'net community productivity';
    case 'nitr'
        FLUXNAME = 'NH_4 nitrification';
    case 'netnitr'
        FLUXNAME = 'net nitrification';
    case 'oxidize'
        FLUXNAME = 'abiotic Fe(II) oxidation';
    case 'pholys'
        FLUXNAME = 'NH_3 photolysis';
    case 'photox'
        FLUXNAME = 'Fe photo-oxidation';
    case 'precip'
        FLUXNAME = 'CaCO_3 precipitation';
    case 'prod'
        FLUXNAME = 'primary production';
    case 'rainout'
        FLUXNAME = 'aerosol rainout';
    case 'reduce'
        FLUXNAME = 'Fe-P reduction';
    case 'revweather'
        FLUXNAME = 'reverse weathering';
    case 'scav'
        FLUXNAME = 'phosphate scavenging';
    case 'sed'
        FLUXNAME = 'sedimentation';
    case 'sfw'
        FLUXNAME = 'seafloor weathering';
    case {'carbsink','sink','bifsink'}
        FLUXNAME = 'sinking';
    case 'sorb'
        FLUXNAME = 'sorption';
    case {'subduct','subd'}
        FLUXNAME = 'subduction';
    case 't'
        FLUXNAME = 'total ocean';
    case 'totremin'
        FLUXNAME = 'total remineralization';
    case 'volc'
        FLUXNAME = 'volcanic outgassing';
    case 'wthr'
        FLUXNAME = 'weathering'; 
    case {'burOC','burIC','burON','burRN'}
        FLUXNAME = 'burial'; 
        SPECIES  = FLUX(4:end); 
    case {'burial','death','export','fixation'}                           % these are already pretty descriptive
        FLUXNAME = FLUX; 
    case {'axrmP','axrmN', 'axrmC','oxrmP','oxrmN','oxrmC','axrm','oxrm'} % total remin, anoxic and oxic
        if strcmp(FLUX(1:4),'axrm')
            FLUXNAME = 'anaerobic remin';
        else
            FLUXNAME = 'aerobic remin';
        end
    case {'SPECIES','species'}                                            % this is a flag that I just want a descriptive species output
        FLUXNAME = ''; 
    otherwise
        error('Input 1 is not a known flux in the model'); 
end
if iscell(desig) % these should be entered in the folowing order: reservoir, species
    switch desig{1}
        case {'a','s','d','n','z','c','u','t'}                            % it do be a reservoir
            if strcmp(desig{1},'t')                                       % replace with 'total'
                desig{1} = 'total';
            end
            RES = desig{1}; SPECIES = desig{2}; 
        otherwise % assume it's a species
            if strcmp(desig{2},'t')                                       % replace with 'total'
                desig{2} = 'total';
            end
            RES = desig{2}; SPECIES = desig{1};
    end
elseif desig == 666                                                       % not inputted
    RES = '';       SPECIES = ''; 
elseif exist('SPECIES','var')                                             % already assigned!
    RES = '';
else % not a cell, decide which it could be
    switch desig
        case {'a','s','d','n','z','c','u','t'}                            % it do be a reservoir
            if strcmp(desig,'t')                                          % replace with 'total'
                desig = 'total';
            end
            RES = desig;    SPECIES = '';
        otherwise
            RES = '';       SPECIES = desig;
    end
end

% give species more clear names
switch SPECIES
    case 'N2'
        SPECIES = 'N_2';
    case 'NH3'
        SPECIES = 'NH_3';
    case 'NH4'
        SPECIES = 'NH_4';
    case {'NO3','HNO3'}
        SPECIES = 'HNO_3';
    case 'H3PO4'
        SPECIES = 'H_3PO_4';
    case 'PO4'
        SPECIES = 'PO_4';
    case {'CaCO3','IC','carb'}
        SPECIES = 'CaCO_3';
    case 'CO2'
        SPECIES = 'CO_2';
    case 'O2'
        SPECIES = 'O_2';
    case 'HCO3'
        SPECIES = 'HCO_3';
    case 'CO3'
        SPECIES = 'CO_3';
    case 'CH4'
        SPECIES = 'CH_4';
    case {'OC','POC','ON','OP'}
        SPECIES = ['Org ',SPECIES(end)];
    case {'org','oxi'}
        SPECIES = 'Org Matter';
    case {'omega','Omega'}
        SPECIES = 'Saturation (\Omega)'; 
    case 'sil'
        SPECIES = 'silicate'; 
    case {'fixN','fixedN','FixedN'}
        SPECIES = 'fixed N'; 
    case 'totN'
        SPECIES = 'N_{total}';
    case 'FeO'
        SPECIES = 'FeO';
    case {'bif','FeOH3'}
        SPECIES = 'Fe(OH)_3'; 
    case 'Fe2O3'
        SPECIES = 'Fe_2O_3'; 
    case {'Feox','ironox'}
        SPECIES = 'reduced Fe'; 
    case {'FeP','FeIIP'}
        SPECIES = 'Fe_3(PO4)_2';
    case 'FePO4'
        SPECIES = 'FePO_4';
    case 'Fe2SiO4'
        SPECIES = 'Fe_2SiO_4';
    case {'C','H','O','N','P','Fe'} % already good!  

end

% and make a descriptive print name out of the two/three inputs
switch order
    case 'rfirst' % reservoir is put first in the name
        fname = [RES,' ',SPECIES,' ',FLUXNAME];                           % ie. == 'd RN diffusion'
    case 'rlast' % reservoir is put last in the name
        fname = [SPECIES,' ',FLUXNAME,' ',RES];                           % ie. == 'RN diffusion d'
    case 'ffirst' % flux is first in the name
        fname = [FLUXNAME,' ',RES,' ',SPECIES];                           % ie. == 'diffusion d RN'
    otherwise 
        error('Input 3 is not a known order'); 
end

end


%% Subfunction : change axis limits
function SetYLim(ax,oyl,lbound,ubound)
    if exist('ubound','var') % if an upper boundary has been pre-defined, apply it 
        set(ax,'ylim',[lbound ubound]); 
    else
        set(ax,'ylim',[lbound oyl(end)]); 
    end
end

