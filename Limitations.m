%% Calculate any limitation using half-saturation uptake parameters
% Julia Horne, 2022

function LIM = Limitations(r,conc,c,v)
    [lim] = DefineProdLimitations(r.s,conc.s,c.s,v); % use subfunction to find limits for primary production
    % define other limitations not limited to surface ocean
    bxs = {'s','d','n','z'};
    for ib = 1:length(bxs)
        bs = bxs{ib}; 
        lim.Hsb.(bs)  = conc.(bs).pH./(conc.(bs).pH + 8);       % phosphate sorption pH sensitivity, [PO4^3-] declines above pH = 7 (Konhauser et al. 2007)
        lim.Mtro.(bs) = c.(bs).O2./(c.(bs).O2 + v.Ks.Omt);      % methanotroph production limited by oxygen level (delta M in Goldblatt et al. 2006)
        lim.Onit.(bs) = c.(bs).O2./(c.(bs).O2 + v.Ks.On);       % nitrification O2 sensitivity (aka V_Onit in EarthN)
        lim.Nnit.(bs) = c.(bs).NH4./(c.(bs).NH4 + v.Ks.Nn);     % nitrification NH4 sensitivity (aka V_NOnit in EarthN)
        lim.Nden.(bs) = c.(bs).HNO3./(c.(bs).HNO3 + v.Ks.Nd);   % denitrification HNO3 sensitivity (aka V_NOde in EarthN)
        lim.Oaer.(bs) = c.(bs).O2./(c.(bs).O2 + v.Ks.Oor);      % aerobic remin O2 sensitivity
        lim.Oden.(bs) = c.(bs).O2./(c.(bs).O2 + v.Ki.Od);       % O2 inhibition of denitrification
        lim.Nana.(bs) = c.(bs).HNO3./(c.(bs).HNO3 + v.Ki.Nm);   % HNO3 inhibition of methanogenesis
        lim.Hnit.(bs) = c.(bs).H./(c.(bs).H + v.Ki.Hn);         % pH inhib nitrification
        lim.Hden.(bs) = c.(bs).H./(c.(bs).H + v.Ki.Hd);         % pH inhib denitrification
        lim.Anox.(bs) = c.(bs).O2./(c.(bs).O2 + v.Ks.Oor);      % degree of anoxicity limitation - defined by functionality of aerobic metabolism
        lim.FeP.(bs)  = c.(bs).FeOH3./(c.(bs).FeOH3 + v.KdP);   % phosphate affinity for Fe-sorption limitation
    end
    
    % check all limitations, and make any that go negative == 0
    lst = fieldnames(lim); 
    for il = 1:length(lst)
       if isstruct(lim.(lst{il}))
          bxs = fieldnames(lim.(lst{il})); 
          for ib = 1:length(bxs)
             LIM.(lst{il}).(bxs{ib}) = CheckNegativeLims(lim.(lst{il}).(bxs{ib})); 
          end
       else 
           LIM.(lst{il}) = CheckNegativeLims(lim.(lst{il})); 
       end
    end
end

%% Subfunction: Calculate primary productivity limitations by nutrients in surface ocean only

function lim = DefineProdLimitations(rS,conS,cS,v) 
    Kn       = 1.6e-6 .* v.conv.L2m3;                           % mol/m^3; half saturation uptake value for NH4 (from 1.6e-6 mol/L; Johnson, 2017)
    Kp       = 0.1e-6 .* v.conv.L2m3;                           % mol/m^3; half saturation uptake value for PO4 (from 0.1e-6 mol/L; Johnson, 2017)
    Knp      = 8;                                               % dimensionless half saturation constant for Redfield sensitivity N2:PO4 (Johnson, 2017)
    Kc       = 4.24e-6 .* v.oc.rho;                             % mol/m3; half saturation uptake of CO2 (3.5-5 µM, Burkhardt et al. 2001)
    
    T0       = 298.15;                                          % K; reference temperature for Henry's Law eqn 
    Tc       = 288;                                             % K; assume a constant temperature of 15C for the purposes of this 'constant'
    p0       = 5000;                                            % Pa; N2 pressure in atm for fixation following Michaelis-Menton eqn   
    K0_N2    = 6.4e-6;                                          % mol/m3Pa; reference Kh at T0 for N2 (Sander, 2015)
    vH_N2    = 1300;                                            % K; van't Hoff parameter for N2 (Sander, 2015)
    Kh_N2    = K0_N2 .* exp(vH_N2.*(1./Tc-1./T0));              % mol/m3Pa; Henry's constant for N2
    KmN2     = Kh_N2*p0;                                        % mol/m3; Michaelis constant (aka K_fix)
    
    totbio_N = rS.RN + rS.HNO3;                                 % mol N; total bioavailable N in surface ocean
    NP_ratio = totbio_N./rS.H3PO4;                              % surface ocean NH3+HNO3/PO4 ratio

    % Calculate the limits according to nutrient availability in surface
    % ocean
    lim.N   = (totbio_N./v.oc.vol.s)./((totbio_N./v.oc.vol.s)+Kn);% dimensionless N limiting factor
    lim.P   = (rS.H3PO4./v.oc.vol.s)./((rS.H3PO4./v.oc.vol.s)+Kp);% dimensionless PO4 limiting factor                            
    lim.C   = (conS.DIC.*v.oc.rho)./((conS.DIC.*v.oc.rho)+Kc);  % Carbon limitation on productivity 
    lim.Phx = (cS.O2./(cS.O2 + v.Ks.Oor));                      % degree of anoxicity limitation on photoferrotrophy -- anoxic zone defined by lower limit of functional aerobic metabolisms (Tiano et al. 2014)
    lim.RN  = cS.RN./(cS.RN + Kn);                              % NH4/NH3 limitation on photoferrotrophy
    lim.RS  = NP_ratio./(NP_ratio + Knp);                       % Redfield sensitivity (proportionates N2 fixation and NH3 assimilation)
    lim.NS  = (rS.N2./v.oc.vol.s)./((rS.N2./v.oc.vol.s)+KmN2);    % atmospheric N2 reservoir sensitivity
end

%% Subfunction: Check that the limitations are not going negative, and make them 0 if they do
function nlim = CheckNegativeLims(LIM)
    nlim = LIM; 
    for it = 1:length(LIM)
        if LIM(it) < 0
            nlim(it) = 0; 
        end
    end
end
