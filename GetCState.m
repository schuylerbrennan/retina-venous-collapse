function [Dc,tauc,Pc,Ac,nc,Lc,muc,M0c,Prmc,Rc,Cc,Qc,...
    condrespsignalc,krogh, rref, etac] = GetCState(Cmeta,Length_constant,multiplier,...
    Pc,tauc,muc,M0c,Prmc,formatvars,respswitches, PIOP_2)


    %Initial guesses for the FEED and MID components(micrometers)
    %%% Initial values of diameter (D) and activation (A) in LA and SA
    Dc.DLAi = 100e-4;
    Dc.DSAi = 45e-4;
    Ac.ALAi = 0.5;
    Ac.ASAi = 0.5;
    
    %%% Time that we will solve the ODEs for
    tvec = [0,300];
    %%% This function contains the ODEs to solve and solves them
    [time,vector] = ode15s(@(t,y) DAnometaDonly(t,y,Pc,Prmc),...
        tvec,[Dc.DLAi Dc.DSAi]);
    
    %  JB-vectorized the for loop using normal matlab matrix manipulations and
    %    a nice function, "deal"
    
    [Dc.DLAArray,Dc.DSAArray] = ...
        deal(vector(:,1)',vector(:,2)');
    
    %Plot diameters vs time
    % figure(100)
    % plot(time, Dc.DLAArray.*10.^4, 'b')
    % xlim([0 50])
    % figure(101)
    % plot(time, Dc.DSAArray.*10.^4, 'r')
    % xlim([0 50])
    
    % returns the last diameters of the solutions
    %%% Get the control state from the solved ODE
    Dc.DLA = Dc.DLAArray(end);
    Dc.DSA = Dc.DSAArray(end);
    
    
    %Invoked symmetry for DLV and DSV in control state
    %%% Assumed capillary diameter and calculated vein diameters from 
    Dc.DC = 6e-4; %assumed capillary diameter
    %%% Calculated from symmetry and eqn 13 in publication 3
    Dc.DSV = (muc.muSV*tauc.tauSA/(muc.muSA*tauc.tauSV))^(1/3)*Dc.DSA;
    Dc.DLV = (muc.muLV*tauc.tauLA/(muc.muLA*tauc.tauLV))^(1/3)*Dc.DLA;
    
    %Flow through a single vessel for LA, SA, caps
    %%% Calculated from Q = (pi/32)*(D^3*tau/mu)
    Qc.QLA = (pi*tauc.tauLA*Dc.DLA^3)/(32*muc.muLA);
    Qc.QSA = (pi*tauc.tauSA*Dc.DSA^3)/(32*muc.muSA);
    Qc.QC = (pi*tauc.tauC*Dc.DC^3)/(32*muc.muC);


        
    %Numbers of vessels based on total calf flow
    %%% Assuming 4 vessels in LA, and use conservation of flow for SA and C
    %%% SV and LV are calculated using symmetry of the vessel
    nc.nLA = 4;
    nc.nSA = nc.nLA*Qc.QLA/Qc.QSA;
    nc.nC = nc.nLA*Qc.QLA/Qc.QC;
    nc.nSV = nc.nSA;
    nc.nLV = nc.nLA; 
    
    %remaining flow
    %%% Calculated using conservation of flow
    Qc.QSV = nc.nLA*Qc.QLA/nc.nSV;
    Qc.QLV = nc.nLA*Qc.QLA/nc.nLV;
    

    Qc.QtotalLA = Qc.QLA*nc.nLA;
    Qc.QtotalSA = Qc.QSA*nc.nSA;

    
    %Conductances
    %%% Calculated from 1/Resistance where
    %%% R = delP/Q
    %%% But we account for total flow, so we multiply by n for that section
    Cc.CLA = nc.nLA*Qc.QLA/Pc.delPla;
    Cc.CSA = nc.nSA*Qc.QSA/Pc.delPsa;
    Cc.CC = nc.nC*Qc.QC/Pc.delPc;
    Cc.CSV = nc.nSV*Qc.QSV/Pc.delPsv;
    Cc.CLV = nc.nLV*Qc.QLV/Pc.delPlv;
    
   
    %Resistances
    %%% 1/Conductance
    Rc.RLA = 1/Cc.CLA;
    Rc.RSA = 1/Cc.CSA;
    Rc.RC = 1/Cc.CC;
    Rc.RSV = 1/Cc.CSV;
    Rc.RLV = 1/Cc.CLV;
    
    %Effecctive lengths
    %%% Poiseuille's Law
    %%% Calculated from L = (pi*R*D^4)/(128*mu)
    %%% This is then converted to not include resistance, so
    %%% L = (pi*D^4*delP)/(128*mu*Q)
    Lc.Leffla = pi*Pc.delPla*Dc.DLA^4/(128*muc.muLA*Qc.QLA);
    Lc.Leffsa = pi*Pc.delPsa*Dc.DSA^4/(128*muc.muSA*Qc.QSA);
    Lc.Leffc = pi*Pc.delPc*Dc.DC^4/(128*muc.muC*Qc.QC);
    Lc.Leffsv = Lc.Leffsa;
    Lc.Lefflv = Lc.Leffla;

    %%% Tension calculated from T = (delP * D)/2 (Law of Laplace)
    Tc.TLA = ((Pc.Pla-Pc.IOP)*Dc.DLA)/2;
    Tc.TSA = ((Pc.Psa-Pc.IOP)*Dc.DSA)/2;
    
    %%% This changes values for pressure, tau, diameter, tension,
    %%% conductance, and resistance to account for an elevated IOP
    


    % M0vector = [1 5 8.2766 10 15 20 30 40 50 60 70];
    % M0vector = [1 8.2766 20];
    
    %  Function call to getscJAcontrol in order to find conducted response
    %  signal (function calculates O2 saturation and ATP levels)
    
   
    %  Info for the plots...quickfix (though not the best) JB
    % save JBploT.Tvars;
    
    %  Currently these are the only outputs used.  Note that the total volume
    %  of the calf is returned in Prmc so it can be used later on when
    %  adding/subtracting vessels.
    
    %%% Calculates and returns control state values for the system running
    [condrespsignalc,krogh,Prmc,systemstatec] = getscJAcontrol(formatvars,Dc,Qc,nc,Lc,M0c,...
        Length_constant,Prmc);
    
    
    %  New coefficients for conducted response (metabolic) term (assume same
    %  for all compartments).  Put these in with other parameters
    
    Prmc.LA.Cmeta = Cmeta;
    Prmc.SA.Cmeta = Cmeta;
    
    
    tissue_CO2_content = CO2converter(Qc.QLA*nc.nLA,systemstatec.endpts.CO2(end));
    tissue_CO2_value = content_to_CO2(tissue_CO2_content);

    %set Young modulus from literature, in dyn/cm^2
    % Prmc.E_vein = 10 * 6e6;
    % Prmc.E_artery = Prmc.E_vein/2;
    % Prmc.E_artery = 2.2e5;
    % Prmc.E_vein =  6.6e5;
    Prmc.E_vein =  3*16e6;
    Prmc.E_artery = 16e6;


    %set h values from literature in cm
    % Prmc.h.LA = 0.00239;
    % Prmc.h.SA = 0.00057;
    % Prmc.h.SV = 0.00034;
    % Prmc.h.LV = 0.00078;

    %Wall to lumen ratio from Guidoboni 2014
    Prmc.h.LA = Dc.DLA*.4;
    Prmc.h.SA = Dc.DSA*.4;
    Prmc.h.SV = Dc.DSV*.05;
    Prmc.h.LV = Dc.DLV*.05;
    Prmc.h.C = Dc.DC*0.17;

    %set nu
    Prmc.nu = 0.49;


    %calculating control state transmural pressure
    Pc.TMLA = Pc.Pla - Pc.IOP;
    Pc.TMSA = Pc.Psa - Pc.IOP;
    Pc.TMC = Pc.Pc - Pc.IOP;
    Pc.TMSV = Pc.Psv - Pc.IOP;
    Pc.TMLV = Pc.Plv - Pc.IOP;

    %solves for rref and etac in the control state
    [rref, etac] = G_rref(Dc, Pc, Prmc);

    %set values for tube laws
    Prmc.kr.LA = 8*pi*muc.muLA;
    Prmc.kr.SA = 8*pi*muc.muSA;
    Prmc.kr.SV = 8*pi*muc.muSV;
    Prmc.kr.LV = 8*pi*muc.muLV;

    Prmc.kp.LA = (Prmc.E_artery)/(12*(1-Prmc.nu^2))*(Prmc.h.LA/rref.LA)^3;
    Prmc.kp.SA = (Prmc.E_artery)/(12*(1-Prmc.nu^2))*(Prmc.h.SA/rref.SA)^3;
    Prmc.kp.SV = (Prmc.E_vein)/(12*(1-Prmc.nu^2))*(Prmc.h.SV/rref.SV)^3;
    Prmc.kp.LV = (Prmc.E_vein)/(12*(1-Prmc.nu^2))*(Prmc.h.LV/rref.LV)^3;

    Prmc.kl.LA = (12*(pi*rref.LA^2))/(pi*Prmc.h.LA^2);
    Prmc.kl.SA = (12*(pi*rref.SA^2))/(pi*Prmc.h.SA^2);
    Prmc.kl.SV = (12*(pi*rref.SV^2))/(pi*Prmc.h.SV^2);
    Prmc.kl.LV = (12*(pi*rref.LV^2))/(pi*Prmc.h.LV^2);

    Prmc.Leffla = Lc.Leffla;
    Prmc.Leffsa = Lc.Leffsa;
    Prmc.Leffc = Lc.Leffc;
    Prmc.Leffsv = Lc.Leffsv;
    Prmc.Lefflv = Lc.Lefflv;

    Prmc.nLA = nc.nLA;
    Prmc.nSA = nc.nSA;
    Prmc.nC = nc.nC;
    Prmc.nSV = nc.nSV;
    Prmc.nLV = nc.nLV;

    

    %% New control state (A = 0.5):

    k.kLA = -respswitches.myo.LA*Prmc.LA.Cmyo*Tc.TLA+...
        respswitches.shear.LA*Prmc.LA.Cshear*tauc.tauLA+...
        respswitches.CO2.LA*Prmc.LA.Cco2*tissue_CO2_value;
    k.kSA = -respswitches.myo.SA*Prmc.SA.Cmyo*Tc.TSA+...
        respswitches.shear.SA*Prmc.SA.Cshear*tauc.tauSA+...
        respswitches.CO2.SA*Prmc.SA.Cco2*tissue_CO2_value;

    %%% Adds metabolic response to Stone equation
    Cpptone.LA = k.kLA+Prmc.LA.Cmeta*condrespsignalc(1);
    Cpptone.SA = k.kSA+Prmc.SA.Cmeta*condrespsignalc(2);
    Prmc.LA.Cpptone = Cpptone.LA;
    Prmc.SA.Cpptone = Cpptone.SA;

    %%% Calculates Activation from Stone equation
    Ac.ALA = 1/(1+exp((-1*respswitches.myo.LA*Prmc.LA.Cmyo*Tc.TLA)+...
        (respswitches.shear.LA*Prmc.LA.Cshear*tauc.tauLA)+...
        Prmc.LA.Cmeta*condrespsignalc(1)+(respswitches.CO2.LA*Prmc.LA.Cco2*tissue_CO2_value)-...
        Prmc.LA.Cpptone));
    Ac.ASA = 1/(1+exp((-1*respswitches.myo.SA*Prmc.SA.Cmyo*Tc.TSA)+...
        (respswitches.shear.SA*Prmc.SA.Cshear*tauc.tauSA)+...
        Prmc.SA.Cmeta*condrespsignalc(2)+(respswitches.CO2.SA*Prmc.SA.Cco2*tissue_CO2_value)-...
        Prmc.SA.Cpptone));

    end 