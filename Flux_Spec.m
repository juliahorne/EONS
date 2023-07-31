%% ============================= EONS Model =============================== 
% Julia Horne, 2018

% Solves for speciation state of carbon and nitrogen dissolved species
% according to current [TA], [DIC], and [H+]. Equations mostly from Zeebe +
% Wolfgladrow (2001).
% The user has the option to assign a different polynomial solver from that
% used in the published EONS model (TBD) - either the original polynomial
% from Z+WG (ZWG), the modified polynomial including NH4 effects on TA (N),
% or using the carbonate_eq7 function from Roberta Hamme (ce7; this needs 
% to be downloaded and included on your path). Default is the modified 
% polynomial - assign a different one as an added 'inp' input 'Poly'

function conc = Flux_Spec(r,T,inp,v)
if ~isfield('Poly','inp')
    Poly = 'N';
else
    Poly = inp.Poly; 
end
spe = {'TA','DIC','RN'}; box = {'s','n','d','z'}; 
Tbx = OceanTemperatures(T); 

for ix = 1:length(box)
    res = box{ix};
    Tbc.(res) = Tbx.(res) - 273.15; % temperaure in celsius
    
    for il = 1:length(spe)
        conc.(res).(spe{il}) = r.(res).(spe{il}) ./ v.oc.m.(res);% mol/kg; concentration
    end

% equilibrium rate constants (Appendix A) 
% K1 : [H+][HCO3-]/[CO2]        = first carbonic acid dissociation constant
% K2 : [H+][CO3^2-]/[HCO3-]     = second carbonic acid dissociation constant
% Kw : [H+][OH-]                = water dissociation constant
% KB : [H+][B(OH)4-]/[B(OH)3]   = boric acid dissociation constant
% Kn : [NH3][H+]/[NH4+]         = ammonium dissociation constant

% NOTATION: 
% K0 = thermodynamic dissolution constant (mol^2/kg^2)
% Kp = pressure-dependant thermodynamic dissolution constant (mol^2/kg^2), 
% K_ = dissociation constants = e^-ln(Kp/K0) == stoichiometric equilibrium/acidity constants (K*)
% k+ = forward reaction rate constant (units vary, either kg/mol s or /s)
% k- = back reaction rate constant (units vary, either kg/mol s or /s)

%% -------------------- Equilibrium Constants -----------------------------
% all equations from appendix A of "CO2 in Seawater: Equilibrium, Kinetics,
% Isotopes" (Zeebe + Wolf-Gladrow, 2001) and are the "recommended"
% formulations using the pH_T (artificial seawater) scale

    % negative ln of first carbonic acid dissociation constant = K01 
        %           H2O + CO2  --K1-->  H+ + HCO3-
    K.K01.(res) = 2.83655 - 2307.1266./Tbx.(res) - 1.5529413.*log(Tbx.(res))...
        - (0.207608410 + 4.0484./Tbx.(res)).*sqrt(v.oc.S.(res))...
        + 0.0846834.*v.oc.S.(res) - 0.00654208.*v.oc.S.(res).^(1.5)...
        + log(1 - 0.001005.*v.oc.S.(res)) ; 
    
    % negative ln of second carbonic acid dissociation constant = K02 
        %             HCO3-  --K2-->  H+ + CO3^2-
    K.K02.(res) = -9.226508 - 3351.6106./Tbx.(res) - 0.2005743.*log(Tbx.(res))...
        - (0.106901773 + 23.9722./Tbx.(res)).*sqrt(v.oc.S.(res))...
        + 0.1130822.*v.oc.S.(res) - 0.00846934.*v.oc.S.(res).^(1.5)...
        + log(1 - 0.001005.*v.oc.S.(res)) ; 
    
    % negative ln of water dissociation constant = K0w 
        %               H2O  --Kw-->  H+ + OH-
    K.K0w.(res) = 148.96502 - 13847.26./Tbx.(res) - 23.6521.*log(Tbx.(res))...
        + (118.67./Tbx.(res) - 5.977 + 1.0495.*log(Tbx.(res))).* sqrt(v.oc.S.(res))...
        - 0.01615.*v.oc.S.(res); 

    % negative ln of Boric acid dissociation constant = K0B 
        %       B(OH)3 + H2O  --Kw-->  B(OH)4- + H+ 
    K.K0B.(res) = (-8966.9 - 2890.53.*v.oc.S.(res).^(0.5)...
        - 77.942.*v.oc.S.(res) + 1.728.*v.oc.S.(res).^(1.5)...
        - 0.0996.*v.oc.S.(res).^2)./Tbx.(res) + 148.0248...
        + 137.1942.*v.oc.S.(res).^(0.5) + 1.62142.*v.oc.S.(res)...
        - (24.4344 + 25.085.*v.oc.S.(res).^(0.5) + 0.2474.*v.oc.S.(res))...
        .*log(Tbx.(res)) + 0.053105.*v.oc.S.(res).^(0.5).*Tbx.(res) ; 

    % negative ln of ammonium dissociation constant = K0n (Millero 95 eq 76)
    K.K0n.(res) = -6285.33./Tbx.(res) + 0.0001635.*Tbx.(res) - 0.25444...
        + (0.46532 - 123.7184./Tbx.(res)).*v.oc.S.(res).^(0.5)...
        + (-0.01992 + 3.17556./Tbx.(res)).*v.oc.S.(res); 

% ** NOTE: equilibrium solubility constants for Calcite/aragonite are in -log10, not -ln! **  
    % negative log of thermodynamic calcite solubility constant = K0c
    K.K0c.(res) = -171.9065 - 0.077993.*Tbx.(res) + 2839.319./Tbx.(res)...
        + 71.595.*log10(Tbx.(res)) + (-0.77712 + 0.0028426.*Tbx.(res)...
        + 178.34./Tbx.(res)).*v.oc.S.(res).^(0.5) - 0.07711.*v.oc.S.(res)...
        + 0.0041249.*v.oc.S.(res).^(1.5); 

    % negative log of thermodynamic aragonite solubility constant = K0a
    K.K0a.(res) = -171.945 - 0.077993.*Tbx.(res) + 2903.293./Tbx.(res)...
        + 71.595.*log10(Tbx.(res)) + (-0.068393 + 0.0017276.*Tbx.(res)...
        + 88.135./Tbx.(res)).*v.oc.S.(res).^(0.5) - 0.10018.*v.oc.S.(res)...
        + 0.0059415.*v.oc.S.(res).^(1.5); 

    % non-pressure sensitive equilibrium constants
    K1np.(res) = exp(K.K01.(res));                          % CO2 acid constant 
    K2np.(res) = exp(K.K02.(res));                          % carbonic acid constant
    Kwnp.(res) = exp(K.K0w.(res));                          % ion water product 
    KBnp.(res) = exp(K.K0B.(res));                          % mol/kg; boric acid constant 
    Knnp.(res) = exp(K.K0n.(res));                          % ammonium acid constant
    Kcnp.(res) = 10.^K.K0c.(res);                           % mol/kg; calcite solubility product 
    Kanp.(res) = 10.^K.K0a.(res);                           % mol/kg; aragonite solubility product

%% check non-pressurized pK values (== -log10(K*)) against surface ocean ests (Zeebe + Wolf-Gladrow, 2001 table A.11.2)
    pK1.(res) = -log10(K1np.(res));                         % == 5.8563 
    pK2.(res) = -log10(K2np.(res));                         % == 8.9249
    pKw.(res) = -log10(Kwnp.(res));                         % == 13.2173
    pKB.(res) = -log10(KBnp.(res));                         % == 8.5975
    pKn.(res) = -log10(Knnp.(res));                         % == 9.19 (from Millero 1995, eq 77)
    pKc.(res) = -log10(Kcnp.(res));                         % == 6.3693
    pKa.(res) = -log10(Kanp.(res));                         % == 6.1883

%% calculate pressure sensitive equilibrium constants a la Z+WG (2001) and their "equic.m" function
    constlist = {'k1','k2','kw','kb','kc','ka'}; Kplist = {'Kp1','Kp2','Kpw','KpB','Kpc','Kpa'}; 
    % unpack pressurization constants for easier calling
    a0 = v.spx.A0.ZWG; a1 = v.spx.A1.ZWG; a2 = v.spx.A2.ZWG; 
    b0 = v.spx.B0.ZWG; b1 = v.spx.B1.ZWG; b2 = v.spx.B2.ZWG; 
    for ik = 1:length(Kplist)
        % A.11 in Zeebe + Wolf-Gladrow, and eqs. 90-92 from Millero 1995
        dV.(Kplist{ik}).(res) = a0.(constlist{ik}) + a1.(constlist{ik}).*Tbc.(res)...
            + a2.(constlist{ik}).*Tbc.(res).^2;             % molal volume change with pressure
        dk.(Kplist{ik}).(res) = b0.(constlist{ik}) + b1.(constlist{ik}).*Tbc.(res)...
            + b2.(constlist{ik}).*Tbc.(res).^2;             % molal compressibility change with pressure
        % negative ln of the pressurized dissociation constants for each
        K.(Kplist{ik}).(res)  = (-dV.(Kplist{ik}).(res) + 0.5 .* dk.(Kplist{ik}).(res) .* v.oc.press.(res))...
            .* (v.oc.press.(res) ./ (v.const.Rhat .* Tbx.(res))); 
    end    
    
%% pressure sensitive equilibrium constants
    K1.(res) = exp(K.Kp1.(res)) .* K1np.(res);              % CO2 acid constant 
    K2.(res) = exp(K.Kp2.(res)) .* K2np.(res);              % carbonic acid constant
    Kw.(res) = exp(K.Kpw.(res)) .* Kwnp.(res);              % ion water product
    KB.(res) = exp(K.KpB.(res)) .* KBnp.(res);              % boric acid constant
    Kc.(res) = exp(K.Kpc.(res)) .* Kcnp.(res);              % calcite solubility constant
    Ka.(res) = exp(K.Kpa.(res)) .* Kanp.(res);              % aragonite solubility constant
    % No given pressure constants for ammonium, so just use the standard
    Kn.(res) = Knnp.(res); 
    
% package equilibrium constants for sending to solvers
    Keq.(res).K1 = K1.(res); Keq.(res).K2 = K2.(res); Keq.(res).Kw = Kw.(res); Keq.(res).KB = KB.(res); Keq.(res).Kn = Kn.(res);     

%% -------------------- Calculate [H] from DIC/Alk ------------------------
switch Poly
    case 'ZWG'          % use original polynomial equation - 5th order
        [h] = PolySolvZWG(conc.(res).DIC,conc.(res).TA,Keq.(res),v.oc.B.(res)); 
        % calculate concentrations in subfunction (mol/kg)
        [conc.(res)] = speciation_concentrations(conc.(res).DIC,conc.(res).TA,conc.(res).RN,h,Keq.(res),v.oc.B.(res));    
      
    case 'N'         % use expanded polynomial, including pH effects from ammonium speciation (TA - NH4) - 6th order
        [h] = PolySolvN(conc.(res).DIC,conc.(res).TA,conc.(res).RN,Keq.(res),v.oc.B.(res)); 
        % calculate concentrations in subfunction (mol/kg)
        [conc.(res)] = speciation_concentrations(conc.(res).DIC,conc.(res).TA,conc.(res).RN,h,Keq.(res),v.oc.B.(res));
      
    case 'ce7'          % run Roberta's carbonate equilibrium code instead (need to have 'carbonate_eq7.m' file on path
        ta = conc.(res).TA .* 1e6; % umol/kg
        dic = conc.(res).DIC .* 1e6; % umol/kg
        [~,ph,co2,hco3,co3,h,oh,boh3,boh4,~,~] = carbonate_eq7(v.oc.S.(res),Tbc.(res),v.oc.press.(res).*10,dic,ta);
        % solve for the rest of the system (mol/kg)
        conc.(res).H = h.*1e-6;    
        conc.(res).CO2 = co2.*1e-6; 
        conc.(res).HCO3 = hco3.*1e-6;
        conc.(res).CO3 = co3.*1e-6;
        conc.(res).OH = oh.*1e-6; 
        conc.(res).BOH4 = boh4.*1e-6; 
        conc.(res).BOH3 = boh3.*1e-6;
        conc.(res).NH3 = conc.(res).RN ./ (1 + conc.(res).H./Kn.(res)); % mol/kg; [NH3]
        conc.(res).NH4 = conc.(res).NH3 .* conc.(res).H ./ Kn.(res); % mol/kg; [NH4]
        conc.(res).pH = ph;

end

%% ------------------------ Carbonate Solubility --------------------------
    % from Zeebe & Wolf-Gladrow (2001) and Millero et al. (1995)
    % modern surface waters Omega_arag = 1-5, should be <1 in deep
    % ocean/seds (like 0.1 to 0.9, not 1e-10!)
    Oc.(res) = (v.oc.Ca.(res) .* conc.(res).CO3) ./ Kc.(res);% calcite saturation state (Omega_calcite)
    Oa.(res) = (v.oc.Ca.(res) .* conc.(res).CO3) ./ Ka.(res);% aragonate saturation state (Omega_aragonite)
    
    % Csat = [CO3] at saturation (Omega = 1) therefore [CO3]sat = Ksp/[Ca]
    Csat.(res) = Ka.(res) ./ v.oc.Ca.(res); 

end

%% ------------------------- Carbon fugacity ------------------------------
% Henry's Law coeff for CO2 (Weiss, 1974) - uses surface temperature rather
% than ocean surface temperature!
KH = exp(-60.2409 + 93.4517 * (100./T) + 23.3585 * log(T./100) + ...
    v.oc.S.s .* (0.023517 - 0.023656 .* (T./100) + 0.0047036 * (T./100).^ 2));% mol/kg 
fCO2 = conc.s.CO2 ./ KH;                            % atm; CO2 fugacity

%% package dissociation constants for output

conc.K1    = K1;       % dissociation constants
conc.K2    = K2; 
conc.Kw    = Kw; 
conc.KB    = KB; 
conc.Kn    = Kn; 
conc.Ka    = Ka;
conc.Kc    = Kc;
conc.KH    = KH; 
conc.pK1   = pK1;       % negative common log of equilibrium constants 
conc.pK2   = pK2; 
conc.pKw   = pKw; 
conc.pKB   = pKB; 
conc.pKn   = pKn; 
conc.pKa   = pKa;
conc.pKc   = pKc; 
conc.Oa    = Oa;       % omega saturation values
conc.Oc    = Oc; 
conc.Csat  = Csat;
conc.fCO2  = fCO2;     % CO2 fugacity 

end


%% Subfunctions: polynomial solvers for H from DIC,TA 

function h = PolySolvZWG(dic,ta,Kq,btot) % the original equation, from Z+WG (2001)
% preallocate vectors
 ply.p0 = zeros(size(dic)); ply.p1 = zeros(size(dic)); ply.p2 = zeros(size(dic)); ply.p3 = zeros(size(dic)); ply.p4 = zeros(size(dic)); ply.p5 = zeros(size(dic));
    % *** NOTE: This equation for H does NOT take into account TA effects from
    % NH4/NH3, so when using this polynomial solver:
        %   [TA] = 2[CO3] + [HCO3] + [OH] + [B(OH)4] - [H]

        % Equation 15, Appendix B (as writen) for surface ocean
        %   DIC(KB + H)(K1*H^2 + 2*K1*K2*H)
        %     = (TA*(KB + H)*H - KB*Btot*H 
        %       - Kw(KB + H) + H^2(KB + H))*(H^2 + K1*H + K1*K2);

        % Create polynomial function components 

        ply.p5 = -1.*ones(size(dic));
        ply.p4 = -ta - Kq.KB - Kq.K1;
        ply.p3 = dic .* Kq.K1 - ta...
            .* (Kq.KB + Kq.K1) + Kq.KB .* btot...
            + Kq.Kw - Kq.KB .* Kq.K1 - Kq.K1.*Kq.K2;
        ply.p2 = dic .* (Kq.KB .* Kq.K1 +...
            2 .* Kq.K1 .* Kq.K2) - ta .* ...
            (Kq.KB .* Kq.K1 + Kq.K1 .* Kq.K2) +...
            Kq.KB .* btot .* Kq.K1 + (Kq.Kw .* ...
            Kq.KB + Kq.Kw .* Kq.K1 - Kq.KB .* ...
            Kq.K1 .* Kq.K2) ;
        ply.p1 = 2 .* dic .* Kq.KB .* Kq.K1 .* Kq.K2...
            - ta .* Kq.KB .* Kq.K1...
            .* Kq.K2 + Kq.KB .* btot .* Kq.K1 ...
            .* Kq.K2 + Kq.KB .* Kq.Kw .* Kq.K1...
            + Kq.Kw .* Kq.K1 .* Kq.K2;
        ply.p0 = Kq.Kw .* Kq.K1 .* Kq.K2 .* Kq.KB;


        [h] = PolySolver(dic,ply); % solve the polynomial

end

function h = PolySolvN(dic,ta,rn,Kq,btot) % expanded polynomial, taking NH4 effects on TA into account
% preallocate vectors
 ply.p0 = zeros(size(dic)); ply.p1 = zeros(size(dic)); ply.p2 = zeros(size(dic)); ply.p3 = zeros(size(dic)); ply.p4 = zeros(size(dic)); ply.p5 = zeros(size(dic)); ply.p6 = zeros(size(dic)); 
%  ** NOTE: all terms are concentrations, except for the speciation
%  constants (Kn, K1, K2, etc) 
% This version of the polynomial subtracts NH4 concentration (derived from
% RN and Kn/H) from TA calculation (NH4 = H*RN/(H + Kn) )
% When using this polynomial solver:
    %   [TA] = 2[CO3] + [HCO3] + [OH] + [B(OH)4] - [H] - [NH4]

%    DIC[ KB*K1*H^3 + 2KB*K1*K2*H^2 + K1*H^4 + 2K1*K2*H^3 + KB*K1*Kn*H^2
%           + 2KB*K1*K2*Kn*H + Kn*K1*H^3 + 2K1*K2*Kn*H^2 ]
%                                   =
%    TA[ KB*Kn*H^3 + Kn*H^4 + KB*H^4 + H^5 + Kn*K1*H^3 + KB*K1*H^3 + K1*H^4
%         + KB*Kn*K1*K2*H + Kn*K1*K2*H^2 + KB*K1*K2*H^2 + K1*K2*H^3 + KB*Kn*K1*H^3 ] 
%    - KB*Bt[ H^4 + K1*H^3 + K1*K2*H^2 + Kn*H^3 + Kn*K1*H^2 + Kn*K1*K2*H ]
%	 + Kw[ KB*H^3 + KB*K1*H^2 + KB*K1*K2*H + H^4 + K1*H^3 + K1*K2*H^2 + Kn*H^3
%         + Kn*K1*H^2 + Kn*K1*K2*H + KB*Kn*H^2 + KB*Kn*K1*H + KB*Kn*K1*K2 ]
%    + KB[ H^5 + K1*H^4 + K1*K2*H^3 ] - H^6 - K1[ H^5 + K2*H^4 ] 
%    - Kn[ H^5 + K1*H^4 + K1*K2*H^3 ] - KnKB[ H^4 + K1*H^3 + K1*K2*H^2 ]
%    - RnKB[ H^4 + K1*H^3 + K1*K2*H^2 ] - Rn[ H^5 + K1*H^4 + K1*K2*H^3]

% left-right yields a 6th order polynomial for [H]
    ply.p6 = -ones(size(dic)); 
    ply.p5 = -(ta + Kq.KB + Kq.K1 + Kq.Kn + rn);
    ply.p4 = dic.*Kq.K1 - ta.*(Kq.KB + Kq.Kn + Kq.K1)...
        + Kq.KB.*btot + Kq.Kw - Kq.KB.*(Kq.K1 + Kq.Kn)...
        - Kq.K1.*(Kq.Kn + Kq.K2) - rn.*(Kq.KB + Kq.K1);
    ply.p3 = dic.*(Kq.K1.*(Kq.KB + Kq.Kn) + 2.*Kq.K1.*Kq.K2)...
        - ta.*(Kq.KB.*Kq.Kn + Kq.KB.*Kq.K1 + Kq.Kn.*Kq.K1 + Kq.K1.*Kq.K2 )...
        + Kq.KB.*btot.*(Kq.Kn + Kq.K1 ) ...
        + Kq.Kw.*(Kq.KB + Kq.Kn + Kq.K1 )...
        - Kq.KB.*Kq.K1.*(Kq.K2 + Kq.Kn )...
        - Kq.Kn.*Kq.K1.*Kq.K2 - rn.*(Kq.K1.*Kq.K2 + Kq.KB.*Kq.K1 ); 
    ply.p2 = dic.*(Kq.K1.*Kq.KB.*Kq.Kn + 2.*Kq.K1.*Kq.K2.*(Kq.KB + Kq.Kn))...
        - ta.*(Kq.KB.*Kq.Kn.*Kq.K1 + Kq.Kn.*Kq.K1.*Kq.K2 + Kq.KB.*Kq.K1.*Kq.K2 )...
        + Kq.KB.*btot.*(Kq.K1.*Kq.K2 + Kq.Kn.*Kq.K1 )...
        + Kq.Kw.*(Kq.KB.*Kq.K1 + Kq.K1.*Kq.K2 + Kq.Kn.*Kq.K1 + Kq.KB.*Kq.Kn )...
        - Kq.KB.*Kq.Kn.*Kq.K1.*Kq.K2 - rn.*Kq.KB.*Kq.K1.*Kq.K2; 
    ply.p1 = 2.*dic.*Kq.KB.*Kq.Kn.*Kq.K1.*Kq.K2...
        - ta.*Kq.KB.*Kq.Kn.*Kq.K1.*Kq.K2...
        + Kq.KB.*btot.*Kq.Kn.*Kq.K1.*Kq.K2...
        + Kq.Kw.*(Kq.KB.*Kq.K1.*Kq.K2 + Kq.Kn.*Kq.K1.*Kq.K2 + Kq.KB.*Kq.Kn.*Kq.K1 ); 
    ply.p0 = Kq.Kw.*Kq.KB.*Kq.Kn.*Kq.K1.*Kq.K2; 

    [h] = PolySolver(dic,ply); % solve the polynomial
end

%% Subfunction: solve chosen polynomial system for [H+]
function h = PolySolver(dic,ply) 
    % solve the polynomial equations at each timestep, rather than trying to
    % feed a giant vector into the polynomial solver (it doesn't appreciate
    % that)
    h = zeros(size(dic)); % initialize vector 
    for ic = 1:length(dic)
        psize = length(fields(ply)); 
        if psize == 6 % use 5th order polynomial solver
            poly = [ply.p5(ic); ply.p4(ic); ply.p3(ic); ply.p2(ic); ply.p1(ic); ply.p0(ic)] ; 
        elseif psize == 7 % use 6th order polynomial solver
            poly = [ply.p6(ic); ply.p5(ic); ply.p4(ic); ply.p3(ic); ply.p2(ic); ply.p1(ic); ply.p0(ic)] ; 
        end
        % use roots function to solve polynomials for [H+]
        roos = roots(poly);
        h(ic) = max(real(roos));
    end
end

%% Subfunction: calculate concentrations from H and Keq 
function cc = speciation_concentrations(dic,ta,rn,h,Kq,btot)
% calculate concentrations of all species!    
    cc.CO2  = dic./(1 + Kq.K1./h + Kq.K1.*Kq.K2./(h.*h));   % mol/kg; [CO2] (eq. 1.1.9)
    cc.HCO3 = dic./(1 + h./Kq.K1 + Kq.K2./h);               % mol/kg; [HCO3] (eq 1.1.10)
    cc.CO3  = dic./(1 + h./Kq.K2 + h.*h./(Kq.K1.*Kq.K2));   % mol/kg; [CO3] (eq 1.1.11) 
    cc.OH   = Kq.Kw ./ h ;                                  % mol/kg; [OH]
    cc.BOH4 = Kq.KB .* btot ./ (Kq.KB + h);                 % mol/kg; [B(OH)4-]
    cc.BOH3 = cc.BOH4 .* h ./ Kq.KB;                        % mol/kg; [B(OH)3]
    cc.NH3  = rn ./ (1 + h./Kq.Kn);                         % mol/kg; [NH3]
    cc.NH4  = cc.NH3 .* h ./ Kq.Kn;                         % mol/kg; [NH4]
    cc.pH   = -log10(h);                                    % pH of reservoir
    % rename output so they aren't overwritten
    cc.TA   = ta; 
    cc.DIC  = dic; 
    cc.RN   = rn; 
    cc.H    = h;  
end

%% Subfunction: calculate ocean box temperatures
function Tbx = OceanTemperatures(Tsurf)
Tbx.s         = Tsurf;                                  % K; vector for surface ocean temperature (equals temp for shallow seds)
Tbx.n         = Tbx.s - 10;                             % K; vector for shallow sed temperature
Tbx.d         = 4+273.15.*ones(size(Tsurf));            % K; vector for temperature of deep ocean water (~4C, most temperature variation occurs in surface ocean, so kept constant)
Tbx.z         = 2.5+273.15.*ones(size(Tsurf));          % K; vector for temperature of deep sediment porewater (~ 2.5C, Tromp et al. 1995)
end

%% ------------- Noted errata in Z+WG (2001) table A.11.1 -----------------
% as quoted below from "equic.m" code provided by Univ.Hawaii
% from lines 480-505 : 
%  " 
%     % index: K1 1, K2 2, Kb 3, Kw 4, Ks 5, Kf 6, Kspc 7, Kspa 8,
%     %        K1P 9, K2P 10, K3P 11
% 
%     %----- note: there is an error in Table 9 of Millero, 1995.
%     %----- The coefficients -b0 and b1
%     %----- have to be multiplied by 1.e-3!
% 
%     %----- there are some more errors! 
%     %----- the signs (+,-) of coefficients in Millero 95 do not
%     %----- agree with Millero 79
% 
%     %----- Millero, 1995 Table 9: KW P-coefficients for water
%     %----- not seawater (cf. Millero, 1983). changed 04/05/14
% 
%     a0 = -[25.5   15.82  29.48  20.02  18.03    9.78  48.76   46. ...
%         14.51 23.12 26.57];
%     a1 =  [0.1271 -0.0219 0.1622 0.1119 0.0466 -0.0090 0.5304  0.5304 ...
%         0.1211 0.1758 0.2020];
%     a2 =  [0.0     0.0   -2.608 -1.409 0.316  -0.942  0.0     0.0 ...
%         -0.321 -2.647 -3.042]*1.e-3;
%     b0 = -[3.08   -1.13   2.84   5.13   4.53    3.91  11.76   11.76 ...
%         2.67 5.15 4.08]*1.e-3;
%     b1 =  [0.0877 -0.1475 0.0    0.0794 0.09    0.054  0.3692  0.3692 ...
%         0.0427 0.09 0.0714]*1.e-3;
%     b2 =  [0.0     0.0    0.0    0.0    0.0     0.0    0.0     0.0 ...
%         0.0 0.0 0.0];
%  "
% NOTE: these errors have been compensated for and the correct values used
% in the ConstantsParameters.m script.