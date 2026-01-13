function [Pnew, taunew, Dnew, Tnew, Cnew, Rnew, etanew, Qnew, kroghnew, condrespsignalnew, Anew, Prmcnew] = Get_Elev_IOP_State(Cmeta, Length_constant, multiplier, Pc, Qc, nc, Lc, Dc, tauc, muc, M0c, ...
    Prmcnew, formatvars, respswitches, PIOP_2, rref)
% 
% Calculates pressures at midpoints and nodes by working backwards from Pv to Pa using Arciero 2013
% Solves for diameters using equation from Arciero 2013 and D=2(rref+η)
% Implements tube laws to solve for resistances by calculating three different constants, Kp, Kr, kL
% Solves for tension using Law of Laplace
% Solves for eta 


  
    %define new capillary diameters, will be the same as GetCState 
    Dnew.DC = 6e-4; %assume this capillary diameter 
    taunew.C = 15;
    
    %set new (will be elevated) IOP, Pa, and Pv
    PIOP_2 = PIOP_2;
    Pnew.IOP = PIOP_2;
    Pnew.Pa = Pc.Pa;
    Pnew.Pv = 15*1333;
    %%%%%Pnew.Pv = 5;
    Pnew.delTot = Pnew.Pa - Pnew.Pv;


    %calculate new flow
    Qnew.QLA = Qc.QLA;
    Qnew.QSA = Qc.QSA;
    Qnew.QC = Qc.QC;
    Qnew.QSV = Qc.QSV;
    Qnew.QLV = Qc.QLV;
    
    %solves for pressure drops, midpoint pressures, and pressures at nodes for elevated
    % IOP, using control from supp Arciero 2013
    Pnew.delPlv = fzero(@(x) (2*rref.LV + (2*(1-Prmcnew.nu^2)*rref.LV^2 * (Pnew.Pv + (x/2) - PIOP_2))/(Prmcnew.E_vein*Prmcnew.h.LV))*x^(1/4)...
        -((128*muc.muLV*Lc.Lefflv*Qnew.QLV)/pi)^(1/4), Pc.delPlv);
    Pnew.Plv = Pnew.Pv + (Pnew.delPlv/2);
    Pnew.Psvlv = Pnew.Pv + Pnew.delPlv;
    
    Pnew.delPsv = fzero(@(x) (2*rref.SV + (2*(1-Prmcnew.nu^2)*rref.SV^2 * (Pnew.Psvlv + (x/2) - PIOP_2))/(Prmcnew.E_vein*Prmcnew.h.SV))*x^(1/4)...
        -((128*muc.muSV*Lc.Leffsv*Qnew.QSV)/pi)^(1/4), Pc.delPsv);
    Pnew.Psv = Pnew.Psvlv + (Pnew.delPsv/2);
    Pnew.Pcsv = Pnew.Psvlv + Pnew.delPsv;
    
    %use old pressure drop in capillaries in new system
    % Pnew.delPc = Pc.delPc;
    % %Pnew.delPc = Pnew.f3*Pnew.delTot;
    % Pnew.Pc = Pnew.Pcsv + (Pnew.delPc/2);
    % Pnew.Psac = Pnew.Pcsv + Pnew.delPc;


    %uses pressure drops percentages from SetCState
    Pnew.f1 = Pc.f1;
    Pnew.f2 = Pc.f2;
    % 
    Pnew.delPla = Pnew.f1*Pnew.delTot;
    Pnew.delPsa = Pnew.f2*Pnew.delTot;
    % 
    % Pnew.Psa = Pnew.Psac + (Pnew.delPsa/2);
    % Pnew.Plasa = Pnew.Psac + Pnew.delPsa;
    % 
    % Pnew.Pla = Pnew.Plasa + (Pnew.delPla/2);
    %Pnew.Panew1 = Pnew.Plasa + Pnew.delPla; %<- check that new Pa matches


    %calculating delPc based off of other dels and total delP
    Pnew.delPc = Pnew.delTot - Pnew.delPla - Pnew.delPsa - Pnew.delPlv - Pnew.delPsv;

    Pnew.Pc = Pnew.Pcsv + (Pnew.delPc/2);
    Pnew.Psac = Pnew.Pcsv + Pnew.delPc;


    Pnew.Psa = Pnew.Psac + (Pnew.delPsa/2);
    Pnew.Plasa = Pnew.Psac + Pnew.delPsa;

    Pnew.Pla = Pnew.Plasa + (Pnew.delPla/2);

    %set new transmural pressures
    Pnew.TM.LA = Pnew.Pla - PIOP_2;
    Pnew.TM.SA = Pnew.Psa - PIOP_2;
    Pnew.TM.C = Pnew.Pc - PIOP_2;
    Pnew.TM.SV = Pnew.Psv - PIOP_2;
    Pnew.TM.LV = Pnew.Plv - PIOP_2;


    
    %new diameters
    Dnew.DLV = 2*rref.LV + (2*(1-Prmcnew.nu^2)*rref.LV^2*Pnew.TM.LV)/(Prmcnew.E_vein*Prmcnew.h.LV);
    Dnew.DSV = 2*rref.SV + (2*(1-Prmcnew.nu^2)*rref.SV^2*Pnew.TM.SV)/(Prmcnew.E_vein*Prmcnew.h.SV);
    Dnew.DLA = 2*rref.LA + (2*(1-Prmcnew.nu^2)*rref.LA^2*Pnew.TM.LA)/(Prmcnew.E_artery*Prmcnew.h.LA);
    Dnew.DSA = 2*rref.SA + (2*(1-Prmcnew.nu^2)*rref.SA^2*Pnew.TM.SA)/(Prmcnew.E_artery*Prmcnew.h.SA);
    
    
    % Tension (Law of Laplace) in vasoactive segment
    Tnew.TLA = (Pnew.Pla-PIOP_2)*Dnew.DLA/2;
    Tnew.TSA = (Pnew.Psa-PIOP_2)*Dnew.DSA/2;
    Tnew.TC = (Pnew.Pc-PIOP_2)*Dnew.DC/2;

% Calculate new flow here??

    %calculating new taus
    taunew.tauLV = 32*muc.muLV*Qnew.QLV/(pi*Dnew.DLV^3);
    taunew.tauSV = 32*muc.muSV*Qnew.QSV/(pi*Dnew.DSV^3);
    taunew.tauSA = 32*muc.muSA*Qnew.QSA/(pi*Dnew.DSA^3);
    taunew.tauLA = 32*muc.muLA*Qnew.QLA/(pi*Dnew.DLA^3);


    %Conductances
    Cnew.CLA = nc.nLA*Qnew.QLA/Pnew.delPla;
    Cnew.CSA = nc.nSA*Qnew.QSA/Pnew.delPsa;
    Cnew.CC = nc.nC*Qnew.QC/Pnew.delPc;
    Cnew.CSV = nc.nSV*Qnew.QSV/Pnew.delPsv;
    Cnew.CLV = nc.nLV*Qnew.QLV/Pnew.delPlv;
   

    %sets constants for tube laws
    Prmcnew.kr.LA = 8*pi*muc.muLA;
    Prmcnew.kr.SA = 8*pi*muc.muSA;
    Prmcnew.kr.SV = 8*pi*muc.muSV;
    Prmcnew.kr.LV = 8*pi*muc.muLV;

    Prmcnew.kp.LA = (Prmcnew.E_artery)/(12*(1-Prmcnew.nu^2))*(Prmcnew.h.LA/rref.LA)^3;
    Prmcnew.kp.SA = (Prmcnew.E_artery)/(12*(1-Prmcnew.nu^2))*(Prmcnew.h.SA/rref.SA)^3;
    Prmcnew.kp.SV = (Prmcnew.E_vein)/(12*(1-Prmcnew.nu^2))*(Prmcnew.h.SV/rref.SV)^3;
    Prmcnew.kp.LV = (Prmcnew.E_vein)/(12*(1-Prmcnew.nu^2))*(Prmcnew.h.LV/rref.LV)^3;


    Prmcnew.kl.LA = (12*(pi*rref.LA^2))/(pi*Prmcnew.h.LA^2);
    Prmcnew.kl.SA = (12*(pi*rref.SA^2))/(pi*Prmcnew.h.SA^2);
    Prmcnew.kl.SV = (12*(pi*rref.SV^2))/(pi*Prmcnew.h.SV^2);
    Prmcnew.kl.LV = (12*(pi*rref.LV^2))/(pi*Prmcnew.h.LV^2);


    %new arteriole  using tube laws
        Rnew.RLA = ((Prmcnew.kr.LA*Lc.Leffla)/(pi*rref.LA^2)^2 * (1 + (Pnew.TM.LA/(Prmcnew.kp.LA*Prmcnew.kl.LA)))^-4)/nc.nLA;
        Rnew.RSA = ((Prmcnew.kr.SA*Lc.Leffsa)/(pi*rref.SA^2)^2 * (1 + (Pnew.TM.SA/(Prmcnew.kp.SA*Prmcnew.kl.SA)))^-4)/nc.nSA;
  
    % %new venules resistances using tube laws
    % if Pnew.TM.SV >= 0
    %     Rnew.RSV = ((Prmcnew.kr.SV*Lc.Leffsv)/(pi*rref.SV^2)^2 * (1 + (Pnew.TM.SV/(Prmcnew.kp.SV*Prmcnew.kl.SV)))^-4)/nc.nSV;
    % else 
    %     Rnew.RSV = ((Prmcnew.kr.SV*Lc.Leffsv)/(pi*rref.SV^2)^2 * (1 - (Pnew.TM.SV/(Prmcnew.kp.SV)))^(4/3))/nc.nSV;
    % end 
    % 
    % if Pnew.TM.LV >= 0
    %     Rnew.RLV = ((Prmcnew.kr.LV*Lc.Lefflv)/(pi*rref.LV^2)^2 * (1 + (Pnew.TM.LV/(Prmcnew.kp.LV*Prmcnew.kl.LV)))^-4)/nc.nLV;
    % else 
    %     Rnew.RLV = ((Prmcnew.kr.LV*Lc.Leffla)/(pi*rref.LV^2)^2 * (1 - (Pnew.TM.LV/(Prmcnew.kp.LV)))^(4/3))/nc.nLV;
    % end 

% %use old tube laws (Poiseuille only)
%     %new venules resistances using tube laws
        Rnew.RSV = ((Prmcnew.kr.SV*Lc.Leffsv)/(pi*rref.SV^2)^2 * (1 + (Pnew.TM.SV/(Prmcnew.kp.SV*Prmcnew.kl.SV)))^-4)/nc.nSV;
        Rnew.RLV = ((Prmcnew.kr.LV*Lc.Lefflv)/(pi*rref.LV^2)^2 * (1 + (Pnew.TM.LV/(Prmcnew.kp.LV*Prmcnew.kl.LV)))^-4)/nc.nLV;

%Note: we are not actually using these, I just left them in from our
%meeting in case you ever wanted to check out these calculations 

    % Ccnew.CLA = nc.nLA*Qnew.QLA/Pnew.delPla;
    % Ccnew.CSA = nc.nSA*Qnew.QSA/Pnew.delPsa;
    % Ccnew.CC = nc.nC*Qnew.QC/Pnew.delPc;
    % Ccnew.CSV = nc.nSV*Qnew.QSV/Pnew.delPsv;
    % Ccnew.CLV = nc.nLV*Qnew.QLV/Pnew.delPlv;
    % 
    % %Resistances
    % %%% 1/Conductance
    % Rcnew.RLA = 1/Ccnew.CLA;
    % Rcnew.RSA = 1/Ccnew.CSA;
    % Rcnew.RC = 1/Ccnew.CC;
    % Rcnew.RSV = 1/Ccnew.CSV;
    % Rcnew.RLV = 1/Ccnew.CLV;



    
    Prmcnew.rref.LA = rref.LA;
    Prmcnew.rref.SA = rref.SA;
    Prmcnew.rref.SV = rref.SV;
    Prmcnew.rref.LV = rref.LV;

    %updates eta based on transmural pressure
    etanew.LV = ((1-Prmcnew.nu^2)*rref.LV^2*Pnew.TM.LV)/(Prmcnew.E_vein*Prmcnew.h.LV);
    etanew.SV = ((1-Prmcnew.nu^2)*rref.SV^2*Pnew.TM.SV)/(Prmcnew.E_vein*Prmcnew.h.SV);
    etanew.LA = ((1-Prmcnew.nu^2)*rref.LA^2*Pnew.TM.LA)/(Prmcnew.E_artery*Prmcnew.h.LA);
    etanew.SA = ((1-Prmcnew.nu^2)*rref.SA^2*Pnew.TM.SA)/(Prmcnew.E_artery*Prmcnew.h.SA);
    
    [condrespsignalnew,kroghnew,Prmcnew,systemstatenew] = getscJAcontrol(formatvars,Dnew,Qnew,nc,Lc,M0c,...
            Length_constant,Prmcnew);
    
    tissue_CO2_contentnew = CO2converter(Qnew.QLA*nc.nLA,systemstatenew.endpts.CO2(end));
    tissue_CO2_valuenew = content_to_CO2(tissue_CO2_contentnew);
    
    knew.kLA = -respswitches.myo.LA*Prmcnew.LA.Cmyo*Tnew.TLA+...
            respswitches.shear.LA*Prmcnew.LA.Cshear*taunew.tauLA+...
            respswitches.CO2.LA*Prmcnew.LA.Cco2*tissue_CO2_valuenew;
    knew.kSA = -respswitches.myo.SA*Prmcnew.SA.Cmyo*Tnew.TSA+...
            respswitches.shear.SA*Prmcnew.SA.Cshear*taunew.tauSA+...
            respswitches.CO2.SA*Prmcnew.SA.Cco2*tissue_CO2_valuenew;

    %%% Adds metabolic response to Stone equation
    Prmcnew.LA.Cpptone = knew.kLA+Prmcnew.LA.Cmeta*condrespsignalnew(1);
    Prmcnew.SA.Cpptone = knew.kSA+Prmcnew.SA.Cmeta*condrespsignalnew(2);

    %%% Calculates Activation from Stone equation
    Anew.ALA = 1/(1+exp((-1*respswitches.myo.LA*Prmcnew.LA.Cmyo*Tnew.TLA)+...
        (respswitches.shear.LA*Prmcnew.LA.Cshear*taunew.tauLA)+...
        Prmcnew.LA.Cmeta*condrespsignalnew(1)+(respswitches.CO2.LA*Prmcnew.LA.Cco2*tissue_CO2_valuenew)-...
        Prmcnew.LA.Cpptone));
    Anew.ASA = 1/(1+exp((-1*respswitches.myo.SA*Prmcnew.SA.Cmyo*Tnew.TSA)+...
        (respswitches.shear.SA*Prmcnew.SA.Cshear*taunew.tauSA)+...
        Prmcnew.SA.Cmeta*condrespsignalnew(2)+(respswitches.CO2.SA*Prmcnew.SA.Cco2*tissue_CO2_valuenew)-...
        Prmcnew.SA.Cpptone));

end