%% ============================= EONS Model =============================== 
% Julia Horne, 2018
% 
% Parameterization of temperature effects of different greenhouse gases.
% optional input PLOT ('on') tells this function to produce the plots throughout;
% do not include this command if just running within the model!
% ****  ----------------  TERMINOLOGY NOTE  ------------------  ****
% y = gamma = optical depth parameters == tau* from McKay et al. (1999)
% (ie. y_ir = thermal optical depth  and y_hz = haze optical depth)

function [tp] = TempParam(v,method,PLOT)
if ~exist('PLOT','var')
    PLOT = 'off'; 
end

%% Temperature parameterization constants calculation
Ts             = 289;              % K; modern surface temp 
pCO2           = v.atm.CO2pal./v.atm.mol; % mixing ratio of CO2 in modern atm (approximately 300 ppm)

%% change in temperature is balanced by with incoming solar flux (Fs),
% outgoing IR flux - blackbody radiation - (Fbb), and the greenhouse
% forcing (GF) from GHGs. 
    % Fbb = sigma Ts^4             Fs = (1-albedo)S/4
    % y_ir = Y(a1 Ts + b1 GF)      Y = (N/N_0)^q         
    % y_hz = a2(CH4/CO2)	       X = (A/M_oc Cp) 
    % dT/dt = X[Fs(1 - y_hz)(1 + 3/4y_ir) + y_hzFs/2 - Fbb]
%  setting dT/dt to 0 makes the blackbody outgoing flux equal to heating
%  fluxes:
    % X Fbb = X[Fs(1 - y_hz)(1 + 3/4y_ir) + y_hzFs/2]
% divide by X (constant) and simplify w/r/to Fs
    % Fbb = Fs[ (1 - y_hz)(1 + 3/4y_ir) + y_hz/2]

%% Define temperature and radiative forcings for Archean and Modern Earth
% for the modern day:
Fs_m = solarrad(v.S.Pref,v);        % W/m2; modern solar flux
Fbb_m = blackbody(Ts,v);            % W/m2; modern black body IR flux
GFm = RFInterp(pCO2,'CO2',v.bb);    % W/m2; modern GHF at pCO2 = 400 ppm

% for the Archean:
Fs_a = solarrad(v.S.Aref,v);        % W/m2; archean solar flux
Fbb_a = Fbb_m;                      % W/m2; archean black body IR flux at modern temp
switch method
    case 'A'
        GFa = RFInterp(v.bb.aPAL*pCO2,'CO2',v.bb); % W/m2; archean GHF at pCO2 (assuming 10^2-10^4 x PAL - Catling + Zanhle, 2020)
    otherwise
        GFa = GFm + 30;             %  W/m2; archean GHF at pCO2 = 15,000 ppm
end

%% Calculate thermal IR opacity (y_ir)
% since we know what the modern day Fs, Fbb, and y_hz = 0 should be, we can
% solve for a reference value y_ir:
    % Fbb = Fs[ (1 - y_hz)(1 + 3/4y_ir) + y_hz/2]
% becomes:
    % Fbb/Fs =  (1 - y_hz)(1 + 3/4y_ir) + y_hz/2
% and eliminate the y_hz term, and assume no pressure broadening effects (Y = 1)
yir_m = ThermalIR(Fbb_m,Fs_m);      % modern thermal IR optical depth

% then we can find constants to equate this to Ts, GF from knowing the
% modern day values for those variables and knowing that Ts(mod) = Ts(arc)
% and GFa
yir_a = ThermalIR(Fbb_a,Fs_a);      % archean thermal IR optical depth

%% Calculate the constants (a1, b1) that equate the two eras
% now solve for the equilibrating constants from:
    % yir_a = a1Ts + b1GF_a     and     yir_m = a1Ts + b1GF_m
    % a1Ts = yir_a - b1GF_a     and     a1Ts = yir_m - b1GF_m 
% since a1Ts is the same for both, equate the two statements:
    % yir_a - b1GF_a = yir_m - b1GF_m
    % b1GF_m - b1GF_a = yir_m - yir_a
% solve for b1:
b1 = bconstant(yir_a,yir_m,GFa,GFm);  % thermal IR optical depth GHG constant

% re-enter into first equation to find a1 for modern Earth (should == Archean)
    % yir_m = a1Ts + b1GFm
a1 = aconstant(yir_m,b1,GFm,Ts,v);% thermal IR optical depth Temp constant

% plot the relationship of yir and greenhouse forcing
switch PLOT
    case 'on'
        gfrange = 0:1:100; % radiative forcings from pCO2
        yir = a1.*Ts + b1.*gfrange; % thermal opacity at this range of GHGs 
        figure(333); clf; 
        plot(gfrange,yir); hold on; 
        ylabel('\gamma_{IR}') ; xlabel('CO_2 RF (W/m^2)'); 
        scatter(GFm,yir_m,'ob','filled'); 
        yira = a1.*Ts + b1.*GFa; % thermal opacity for archean GHG
        scatter(GFa,yira,'or','filled'); 
        legend('\gamma_{IR} according to changing GF so T_s = 289 K','GF = 35 W/m^2',['GF = ',num2str(GFa,'%.f'),' W/m^2'],...
            'location','southeast','NumColumns',2); 
end

%% calculate q 
% % From Goldblatt et al. (2009) : doubling PAN @ 2.5 Ga = +4.4C
% Y = (N/N_0)^q  --> the portion of incoming solar flux retained by N-pressure broadening of absorption bands
% with no haze effects:
    % Fbb/Fs =  1 + 3/4 y_ir
% where:
    % y_ir = Y(a1 Ts + b1 GF)
% therefore:
    % Fbb/Fs =  1 + 3/4[ (N/N_0)^q (a1 Ts + b1 GF) ]
% using this relationship with y_ir, calculate constant q by assuming Ts + 4.4C under 2PAN 
% then we can solve for q:
    % yir2pan = (2N/N_0)^q (a1Ts + b1GF)
    % yir2pan/(a1Ts + b1GF) = (2N/N_0)^q
    % yir2pan/(a1Ts + b1GF) = 2^q
% take the log of both sides:
    % log(yir2pan/(a1Ts + b1GF) ) = qlog(2)

% calcualte the Y and q values at 2PAN in the Archean (ie. Gamma2PAN)
[qa,Gamma2PANa] = n2press(Ts,Fs_a,a1,b1,GFa,v);
% do the same for a modern solar constant and GHG level
[qm,Gamma2PANm] = n2press(Ts,Fs_m,a1,b1,GFm,v);


% plot Gamma w/r/to pN2 for Archean, modern eras
switch PLOT
    case 'on'
        rN2 = logspace(18,21,20);           % N2 atm reservoir
        fN2 = rN2./v.atm.N2pal;             % N2 PAL
        qval = [qa, qm];                    % archean and modern calculated q exponents
        eranms = {'archean','modern'}; refgammas = [Gamma2PANa, Gamma2PANm];
        figure(); clf; hold on; grid on; 
        for iv = 1:2
            GammaN{iv} = fN2.^qval(iv);     % pressure broadening factor for either era
            plot(fN2,GammaN{iv},'DisplayName',['\Gamma_N ',eranms{iv}]); 
            scatter(2,refgammas(iv),'or','filled','DisplayName',[eranms{iv},' 2PAN \Gamma']);
        end
        set(gca,'xscale','log');
        xlabel('f_{N2} (mol N2/mol N2_{ref})');  ylabel('\Gamma_N = {f_{N2}}^q'); 
        legend('-DynamicLegend','location','northwest'); 
        title('N_2 pressure broadening parameterizations');  

        % recreate temperature w/r/to pCO2 at various N levels (fig 1 from Goldblatt et al. 2009)
        pco2 = logspace(-4,-1,10);      % bar; pCO2 
        n2 = [0.5.*v.atm.N2pal, v.atm.N2pal, 2.*v.atm.N2pal,3.*v.atm.N2pal];% mol; atm N2
        inp.era = 'A' ; inp.solcon = v.S.Aref; 
        symb = {'-.','--','-',':'};
        figure(); clf; inp.component = 'dumb'; 
        for in = 1:length(n2) 
            N(in) = (n2(in)./v.atm.N2pal).^qa;   % pressure broadening effect from pN2 in Archean
            yirx{in} = []; Tsurf{in} = []; 
            for ip = 1:length(pco2)
                co2(ip) = pco2(ip).*v.atm.mol;  % mol; atm CO2
                [rf(ip),~,Tc(ip)] = Flux_TempEq(0,0,co2(ip),CH4_0.*v.atm.mol,inp,v);
                % calculate thermal IR opacity (y_ir) with pressure broadening parameter
                yirx{in}(ip) = N(in).*(a1.*Tc(ip).t + (rf(ip).CO2+rf(ip).CH4).*b1);
                Tsurf{in}(ip) = calctemp(Fs_a,yirx{in}(ip),0,v); 
            end
            plot(pco2,Tsurf{in}-273,symb{in},'DisplayName',['fN2=',num2str(n2(in)/v.atm.N2pal)]); hold on; 
        end
        set(gca,'xscale','log');grid on;
        scatter(1.5e-2,15,'ro','filled','DisplayName','Archean pCO2 for Modern Temp'); 
        xlabel('pCO2 (bar)'); ylabel('T_{surf} (^{\circ}C)'); 
        legend('-DynamicLegend','location','northwest','NumColumns',2); 
        title('Pressure Broadening effect at different pCO2 levels in Archean');
end

%% calculate y_hz constants
% From McKay et al. (1999) and DeWitt et al. (2009): y_hz = 1-exp(-B g y_ir)
% where 0.28 < g < 0.5 (for Earth-like hazes, conservative est = 0.5) 
% and B = 1.5 (average enhancement factor from Trainer et al. 2006)
% optical depth of 0.55 = antigreenhouse effect = y_hz when CH4:CO2 > 0.1

% calculate a range of haze factors (y_hz) from a range of thermal
% opacities (y_ir == gamma_ir)
gamma_ir = 0:1:10;                  % infrared opacity
gamma_hz = hazefactor(gamma_ir,v);  % antigreenhouse parameter from 0-10 y_ir
MKT = calctemp(Fs_m,gamma_ir,gamma_hz,v); % surface temperature from McKay et al. 1999

% Recreate DeWitt et al. 2009 fig. 6
switch PLOT
    case 'on'
        figure(666); clf; 
        plot(gamma_ir,MKT,'DisplayName','T_s given g = 0.5 with haze');hold on; 

        % using the equation for the antigreenhouse param from McKay et al. 1999
        % doesn't require finding constants a2, b2 for y_hz:
        yhz_a =  hazefactor(yir_a,v);       % fraction reflected by antigreenhouse haze Archean
        yhz_m =  hazefactor(yir_m,v);       % fraction reflected by antigreenhouse haze modern

        % let's see what the temperature looks like with these calculations
        TEMPHAZY = calctemp(Fs_m,yir_m,yhz_m,v); 
        TEMPNOHZ = calctemp(Fs_m,yir_m,0,v); 

        scatter(yir_m,TEMPHAZY,'ob','filled','DisplayName',...
           sprintf(['T_s given g = 0.5, y_{IR} = ',num2str(yir_m,'%4.2f'),' with haze']));
        scatter(yir_m,TEMPNOHZ,'or','filled','DisplayName',...
            sprintf(['T_s given g = 0.5, y_{IR} = ',num2str(yir_m,'%4.2f'),' without haze']));
        xlabel('\gamma_{IR} (infrared opacity)'); 
        ylabel('T_{surf} (K)'); 

        % tune to a proper value for "g", since 0.5 is a lil high
        glist = 0.2:0.05:0.5;                      % values for g (equal to gamma_ir size)
        for ig = 1:length(glist)
            yhz_vg = 1-exp(-v.atm.hz.beta.*glist(ig).*gamma_ir); % y_hz for g value
            temp_vg = calctemp(Fs_m,gamma_ir,yhz_vg,v); 
            plot(gamma_ir,temp_vg,':','DisplayName',sprintf(['T_s given g = ',num2str(glist(ig)),' with haze'])); 
        end

        % show some other values for y_hz at my calculated y_ir value, from
        % different values of g, and note the change from a no haze temperature
        % with this parameterization
        gval = [0.28, 0.3, 0.35]; sym = {'dr','*r','pr'};
        for igv = 1:length(gval)
            yg(igv)  = 1-exp(-v.atm.hz.beta.*gval(igv).*yir_m);   % haze factor for g value chosen
            ygT(igv) = calctemp(Fs_m,yir_m,yg(igv),v);            % temp for g value chosen
            dgT(igv) = TEMPNOHZ - ygT(igv);                       % difference between no haze and haze at this parameterization
            scatter(yir_m,ygT(igv),sym{igv},'DisplayName',sprintf(['T_s given g = ',num2str(gval(igv)),' y_{IR} = ',num2str(yir_m,'%4.2f'),' with haze, dT = ',num2str(dgT(igv),'%.f'),' K']));

        end
        legend('-DynamicLegend','location','southoutside','NumColumns',2); set(gca,'xlim',[0 3]); 
        title({'DeWitt et al. (2009) Fig. 6';'using McKay et al. (1999) temperature function'}); 
        printplotpdf(20,20,'Tparam_GammaIRvsT'); 
end

% we use the same assumptions as Dewitt (2009) and Trainer (2006) for
% behavior of haze on early Earth - the conservative limit of g = 0.5 (only
% moderately reflective haze) used by Dewitt et al, the scaling factor for 
% haze at CH4:CO2 = 1 used by Dewitt (2009) and Trainer (2006; beta), and
% Dewitt's formation. 
yhz = 1-exp(-1.5.*0.5.*yir_a);   % antigreenhouse factor on early Earth
% therefore:
    % a2 = y_hz./(CH4/CO2) and Dewitt uses beta = 1.5 at CH4:CO2 = 1, so:
a2 = yhz; 

%% export all new constants within a structure
tp.aT  = a1; 
tp.ahz = a2; 
tp.bT  = b1;
tp.q   = qa; % the Goldblatt et al. 2009 paper was studying the Archean, so use this reference value! 
tp.CS  = climatesensitivity(pCO2,tp,v);    % find climate sensitivity (dK per CO2 doubling, expect 2.5-4)

end

%% subfunctions for different RF and Temp equations

function Fbb = blackbody(T,v)
% planets lose thermal radiation from the surface at a level dependent on
% the surface temperature:
    % Fbb = sigma Ts^4 
Fbb = v.const.bol .* T.^4;                          % W/m2, outgoing flux

end

function Fsol = solarrad(S,v)
% incoming solar radiation is limited by planetary albedo and the portion
% of the surface that the radiation hits (~1/4 total surface area at any time)
    % Fs = (1-albedo)S/4
Fsol = (1 - v.ea.alb) .* S .* (1/4);                % W/m2, incident flux
end

function yir = ThermalIR(Fbb,Fs)
% black body radiation (Fbb, outflux) and solar flux (Fs, influx) relate to 
% change in surface temperature, depending on the thermal capacity of the
% atmosphere (yir) and ocean heat capacity (X) in a *simple* in-out function:
    % dT/dt = X[Fs( (1 - y_hz)(1 + 3/4y_ir) + y_hz/2) - Fbb]
%  setting dT/dt to 0 makes the blackbody outgoing flux equal to heating
%  fluxes:
    % X Fbb = X [Fs( (1 - y_hz)(1 + 3/4y_ir) + y_hzFs/2)]
% assuming a temperature change dT/dt = 0, and assuming no haze (y_hz == 0) 
% we solve for yir:
    % Fbb = Fs ( (1 - 0)(1 + 3/4y_ir) + 0/2)
% simplified temperature change without haze
    % Fbb = Fs (1 + 3/4y_ir) 
    % y_ir = 4/3 (( Fbb / Fs) - 1)
    yir = (4/3) .* ((Fbb./Fs) - 1); % fraction (IR opacity)
end

function acon = aconstant(yir,bcon,RF,T,v)
% simplified to exclude pressure broadening (Y == gamma)
    % y_ir = Y(a1 Ts + b1 RF) = a1 Ts + b1 RF
    % a1 = (y_ir - b1 RF) / Ts
    
    acon = (yir - (bcon.*RF))./ (T - v.atm.Tx); % relation constant with supergreenhouse effect
end

function bcon = bconstant(yirA,yirM,RFA,RFM)
% solve for the equilibrating constants between archean (a) and modern (m):
    % yir_a = a1Ts + b1GF_a     and     yir_m = a1Ts + b1GF_m
    % a1Ts = yir_a - b1GF_a     and     a1Ts = yir_m - b1GF_m 
% since a1Ts is the same for both (ie. same surface temp at different RF), equate the two statements:
    % yir_a - b1GF_a = yir_m - b1GF_m
    % b1GF_m - b1GF_a = yir_m - yir_a
    bcon = (yirM - yirA) ./ (RFM - RFA); % relation constant
end

function [q,gamma] = n2press(T,Fs,acon,bcon,RF,v)
% From Goldblatt et al. (2009) : doubling PAN @ 2.5 Ga = +4.4C
T2pan = T + 4.4; % K
% Y = (N/N_0)^q  --> the portion of incoming solar flux retained by
% N-pressure broadening of absorption bands (gamma = Y)
% with no haze effects causes the same temp:
    % Fbb/Fs =  1 + 3/4 y_ir
% pressure broadening modifies thermal opacity, such that:
    % y_ir = Y(a1 Ts + b1 GF)
% therefore:
    % Fbb/Fs =  1 + 3/4[ (N/N_0)^q (a1 Ts + b1 GF) ]
    
%% first, calculate black body radiation under 2 PAN 
Fbb2pan = v.const.bol .* (T + 4.4).^4; % W/m2
% and calculate the thermal opacity under this atm pressure (yir2pan)
yir2pan = ThermalIR(Fbb2pan,Fs);

%% next, calculate constant q (exponent)
% % using this relationship with y_ir, we can solve for q:
%     % yir2pan = (2N/N_0)^q (a1Ts + b1GF)
%     % yir2pan/(a1Ts + b1GF) = (2N/N_0)^q
% % since we only double N2 pressure, then 2N/N_0 == 2:
%     % yir2pan/(a1Ts + b1GF) = 2^q
% % take the log of both sides:
%     % log(yir2pan/(a1Ts + b1GF) ) = qlog(2) 
% % therefore:
%     % q = log(yir2pan / a1Ts + b1GF) / log(2)
%   
Tv = (T2pan - v.atm.Tx); 
q = log(yir2pan ./ (acon.*Tv + bcon.*RF))./log(2); 


%% OR: assume the 4.4 increase is true for both eras, and the q factor equilibrates this effect
% % then dT/dt = 4.4 for both Archean and Modern when at 2PAN, therefore:
%     % dT/dt_m = X[Fs_m (1 + 3/4 y_irm)  - Fbb] == 4.4 == dT/dt_a = X[Fs_a (1 + 3/4 y_ira)  - Fbb]
%     % X[Fs_m (1 + 3/4 y_irm)  - Fbb] = X[Fs_a (1 + 3/4 y_ira)  - Fbb] 
% % assuming same surface temperature, Fbb is equivalent. We assume X is also
% % equivalent, giving:
%     % Fs_m (1 + 3/4 y_irm) =  Fs_a (1 + 3/4 y_ira) 
%     % Fs_m y_irm = Fs_a y_ira
% % and the effects of pressure broadening for 2PAN are simplified to:
%     % Y = 2^q
% % applying that to the thermal IR opacity:
%     % y_ir = Y(a1*Ts + b1*RF) == 2^q (a1*Ts + b1*RF)
% % So we can make the equivalence between Archean broadening and Modern, and solve for q:
%     % Fs_m [2^q (a1*Ts + b1*RFm)] = Fs_a [2^q (a1*Ts + b1*RFa)]
%     % 2^q = Fs_a (a1*Ts + b1*RFa) / Fs_m (a1*Ts + b1*RFm)
% 
% q = log(FsA.*(acon.*T + bcon.*RFA) ./ FsM .*(acon.*T + bcon.*RFM)) ./ log(2) ;


%% finally, calculate Y (gamma) for this doubled N2 pressure
    % Y = (N/N_0)^q 
gamma = 2.^q; % constant

end

function yhz = hazefactor(yir,v)
% From McKay et al. (1999) and DeWitt et al. (2009): y_hz = 1-exp(-B g y_ir)
% where 0.28 < g < 0.5 (for Earth-like hazes, conservative est = 0.5) 
% and B = 1.5 (average enhancement factor from Trainer et al. 2006)
% optical depth of 0.55 = antigreenhouse effect = y_hz when CH4:CO2 > 0.1

hz.lim     = 1e14.*1e-3./v.const.mm.C;                  % mol C/yr; organic haze production minimum for shielding (Wolf + Toon, 2010) ~ 8.3e12
hz.aGH     = 6.8108e14;                                 % mol C/yr; organic haze production threshold for generating anti-greenhouse effect (using Wolf + Toon, 2010 equation 8 with CH4:CO2 ~ 0.2) 
hz.beta    = 1.5;                                       % haze production average enhancement parameter (Trainer et al. 2006)
hz.gref    = 0.5;                                       % haze visibility thickness 0.28 < g < 0.5 (for Earth-like hazes, conservative est = 0.5) DeWitt et al 2009 and McKay et al. 1999
hz.g       = 0.28;                                      % haze visibility thickness tuned to make it so that hazes don't decrease surface temp more than 15K
yhz        = 1 - exp(-hz.beta .* hz.gref .* yir);       % antigreenhouse parameter

end

function dTdt = tempchange(T,S,yhz,PAN,a,b,RF,v)
% temperature change is the balance between incoming (Fs, solar influx) and 
% retained (yir, thermal opacity, for greenhouse gases) radiation, an 
% antigreenhouse factor (yhz, IR reflective capacity), and outgoing thermal
% flux (blackbody radiation, Fbb). Also depends on the heat capacity of the
% ocean (X). 
    % dT/dt = X[Fs ((1 - y_hz)(1 + 3/4y_ir) + y_hz/2) - Fbb]

% The user has the option to pre-calculate a greenhouse factor, or turn it
% off completely. This also is true for pressure broadening, which is
% determined if PAN > 1 (if PAN == 1, then nyir == yir). 

%% first, calculate all components of the equation, given the inputs
X     = v.ea.sa./(v.oc.m.t.*v.oc.hcap);     % Km2/Wyr; total ocean heat capacity per year
Fbb   = blackbody(T,v);                     % W/m2; outgoing IR flux, based on temperature
Fsol  = solarrad(S,v);                      % W/m2; incoming IR flux, based on solar constant
yir   = thermalIR(Fbb,Fsol);                % thermal IR opacity, based on incoming-outgoing flux, with NO pressure broadening
[q,~] = n2press(T,Fsol,a,b,RF,v);           % N2 pressure broadening factor, based on radiative forcing and solar flux
nyir  = (PAN.^q).*(a.*T + b.*RF);           % thermal IR opacity, given pressure broadening factor (gamma = N/N_0^q)

%% now calculate temperature change
dTdt  = X .* (Fsol .* ( (1-yhz) .* (1 + (3/4).*nyir) + (yhz./2) ) - Fbb); % K/yr; temperature change

end

function TMP = calctemp(Fs,yir,yhz,v)
% temperature at different haze, solar influx, and radiative forcing levels
% is found from the Stephan-Boltzmann law, neglecting any geothermal
% heating term (Fg):
    % sigma T^4 = Fs( (1 - y_hz)(1 + 3/4y_ir) + y_hz/2 ) 
    % T = [ Fs( (1 - y_hz)(1 + 3/4y_ir) + y_hz/2 ) / sigma]^(1/4)

TMP = ( (Fs .* ((1-yhz).*(1 + (3/4).*yir) + (yhz./2))) ./ v.const.bol).^(1/4) ;

end

function CS = climatesensitivity(CO2,tp,v)
% With our prototype parametrization, find the climate sensitivity by doubling CO2

GF_0 = RFInterp(CO2,'CO2',v.bb) ;        % greenhouse forcing modern CO2
GF_d = RFInterp(2.*CO2,'CO2',v.bb) ;     % greenhouse forcing with doubled CO2
gfs = {GF_0, GF_d};

for ig =1:length(gfs)
   Teq{ig} = equilibtemp(gfs{ig},v);     % find equilib temp
   yir{ig} = (tp.bT .* gfs{ig}) + (tp.aT .* (Teq{ig} - v.atm.Tx) ); % thermal ir
   nT{ig} = calctemp((v.S.Pref./4).*(1-v.ea.alb),yir{ig},0,v); % new temperature!
end

function eqT = equilibtemp(RF,v)
    eqT = (((((v.S.Pref/4).*(1-v.ea.alb))+(v.bb.a.*RF.^v.bb.b))./v.const.bol).^(1/4)); % equilibrium temp at this RF from pCO2
end

CS = nT{2} - nT{1};                      % high minus low pCO2 temp change

end
