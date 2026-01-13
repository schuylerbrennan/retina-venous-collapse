function [Pvar,Qvar,tauvar,Rvar,Cvar] =  obtainnewflowquantselevtt(Dvar,Dc,Pc,nvar,...
    Lc,muc,Rc,Cc,Qc, Prmc, ind, rref, IOPin, Pvar)


%  Given new diameters and old flow information from the previous control
%  state, this program finds the new flow information
%
%  Inputs--all are structures with fields as described below
%    Dvar-new diameters that we want the flow info for
%    Pc-Previous pressures, we use this to obtain the total pressure drop
%      in the system
%    Rc-Previous resistances...only three of these change
%    nvar-number of vessels in each compartment (may or may not vary
%      depending on version of code)
%    muc-the flow resistances...this may eventually change
%

% Calculates new flows and resistances
% Uses tube laws and transmural pressure to calculate resistance, but veins have fixed diameters 
% 
% Pvar = Pc;
Pvar.Pv = 15*1333;
Pvar.delTot = Pvar.Pa - Pvar.Pv;



% Calulates physical constants
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

Prmc.kr.C = 8*pi*muc.muC;
Prmc.kp.C = (Prmc.E_vein)/(12*(1-Prmc.nu^2))*(Prmc.h.C/rref.C)^3;
Prmc.kl.C = (12*(pi*rref.C^2))/(pi*Prmc.h.C^2);



%  Initialize Rvar and Cvar structure
Rvar = Rc;
Cvar = Rc;

% Solves for radial deformation
etavar.LA = (Dvar.DLA/2) - rref.LA;
etavar.SA = (Dvar.DSA/2) - rref.SA;


%Resistances using Poiseuille'e Law
Rvar.RSA = (128*Lc.Leffsa*muc.muSA)/(pi*nvar.nSA*Dvar.DSA^4);
Rvar.RLA = (128*Lc.Leffla*muc.muLA)/(pi*nvar.nLA*Dvar.DLA^4);
Rvar.RC = (128*Lc.Leffc*muc.muC)/(pi*nvar.nC*Dvar.DC^4);

% solves for arteriole alpha (A/A_ref)
Rvar.alphaLA = (pi*(Dvar.DLA/2)^2)/(pi*rref.LA^2);
Rvar.alphaSA = (pi*(Dvar.DSA/2)^2)/(pi*rref.SA^2);

% Sets tolerance for resistance convergence loop
Rtol = 1e-6;
mydiff = 2*Rtol;

while mydiff > Rtol
    
    % Sets Rvarprev to current resistance values
    Rvarprev = Rvar;
    
    % Conductance calculations
    Cvar.CLA = 1/Rvar.RLA;
    Cvar.CSA = 1/Rvar.RSA;
    Cvar.CC = 1/Rvar.RC;
    Cvar.CSV = 1/Rvar.RSV;
    Cvar.CLV = 1/Rvar.RLV;
    
    
    % %Flows
    % %%% Calculate total Resistance by taking sum of 1/Conductance
    Qvar.QLAden = 1/Cvar.CLA + 1/Cvar.CSA + 1/Cvar.CC + 1/Cvar.CSV + 1/Cvar.CLV; %sum of resistances
    %%% Calculate total flow by taking pressure drop over resistance
    Qvar.Qtotal = (Pvar.Pa-Pvar.Pv)/Qvar.QLAden;
    % %%% Calculate flow for each section of vessel by dividing total flow by
    % %%% number of vessels in section
    Qvar.QLA = Qvar.Qtotal/nvar.nLA;
    Qvar.QSA = Qvar.Qtotal/nvar.nSA;
    Qvar.QC = Qvar.Qtotal/nvar.nC;
    Qvar.QSV = Qvar.Qtotal/nvar.nSV;
    Qvar.QLV = Qvar.Qtotal/nvar.nLV;
    % 
    
    %%% Calculate Pressure at each node as well as midpoint pressure using
    %%% Ohm's Law
    Pvar.Plasa = Pvar.Pv + Qvar.Qtotal*(1/Cvar.CLV + 1/Cvar.CSV + 1/Cvar.CC + 1/Cvar.CSA);
    Pvar.PLA = 0.5*(Pvar.Pa + Pvar.Plasa);
    Pvar.Psac = Pvar.Pv + (Qvar.Qtotal*(1/Cvar.CLV +  1/Cvar.CSV + 1/Cvar.CC));
    Pvar.PSA = 0.5*(Pvar.Plasa + Pvar.Psac);
    Pvar.Pcsv = Pvar.Pv + Qvar.Qtotal*(1/Cvar.CLV + 1/Cvar.CSV);
    Pvar.PC = 0.5*(Pvar.Psac + Pvar.Pcsv);
    Pvar.Psvlv = Pvar.Pv + Qvar.Qtotal*(1/Cvar.CLV);
    Pvar.PSV = 0.5*(Pvar.Pcsv + Pvar.Psvlv);
    Pvar.PLV = 0.5*(Pvar.Psvlv + Pvar.Pv);

    % Calculates compartment pressure drops using node pressures
    Pvar.delPlv = Pvar.Psvlv - Pvar.Pv;
    Pvar.delPsv = Pvar.Pcsv - Pvar.Psvlv;
    Pvar.delPc = Pvar.Psac - Pvar.Pcsv;
    Pvar.delPsa = Pvar.Plasa - Pvar.Psac;
    Pvar.delPla = Pvar.Pa - Pvar.Plasa;

    % Calculates transmural pressures
    Pvar.TMLA = Pvar.PLA - IOPin;
    Pvar.TMSA = Pvar.PSA - IOPin;
    Pvar.TMC = Pvar.PC -  IOPin;
    Pvar.TMSV = Pvar.PSV -  IOPin;
    Pvar.TMLV = Pvar.PLV -  IOPin;

    % Calculates percent of the total pressure drop that occurs in each
    % compartment
    Pvar.f1 = (Pvar.Pa-Pvar.Plasa)/Pvar.delTot;
    Pvar.f2 = (Pvar.Plasa-Pvar.Psac)/Pvar.delTot;
    Pvar.f3 = (Pvar.Psac-Pvar.Pcsv)/Pvar.delTot;
    Pvar.f4 = (Pvar.Pcsv-Pvar.Psvlv)/Pvar.delTot;
    Pvar.f5 = (Pvar.Psvlv-Pvar.Pv)/Pvar.delTot;

    % Calculates venule diameters using equations from Cassani thesis (not
    % currently used)
    Dvar.DLV = 2*rref.LV + (2*(1-Prmc.nu^2)*rref.LV^2*Pvar.TMLV)/(Prmc.E_vein*Prmc.h.LV);
    Dvar.DSV = 2*rref.SV + (2*(1-Prmc.nu^2)*rref.SV^2*Pvar.TMSV)/(Prmc.E_vein*Prmc.h.SV);
    % Dvar.DC = 2*rref.C + (2*(1-Prmc.nu^2)*rref.C^2*Pvar.TMC)/(Prmc.E_vein*Prmc.h.C);
    % Dvar.DLA = 2*rref.LA + (2*(1-Prmc.nu^2)*rref.LA^2*Pvar.TMLA)/(Prmc.E_artery*Prmc.h.LA);
    % Dvar.DSA = 2*rref.SA + (2*(1-Prmc.nu^2)*rref.SA^2*Pvar.TMSA)/(Prmc.E_artery*Prmc.h.SA);
    
    % Updates tau calculations
    tauvar.tauLA = 32*muc.muLA*Qvar.QLA/(pi*Dvar.DLA^3);
    tauvar.tauSA = 32*muc.muSA*Qvar.QSA/(pi*Dvar.DSA^3);
    tauvar.tauC = 32*muc.muC*Qvar.QC/(pi*Dvar.DC^3);
    
    
      % new tube laws for venous resistance
    %   if Pvar.TMSV >= 0
    %       Rvar.RSV = (((Prmc.kr.SV*Lc.Leffsv)/((pi*rref.SV^2))^2)*(1+(Pvar.TMSV/(Prmc.kp.SV*Prmc.kl.SV)))^-4)/nvar.nSV;
    %       Rvar.alphaSV = ((Pvar.TMSV)/(Prmc.kp.SV*Prmc.kl.SV) + 1)^2;
    %       Rvar.ASV = Rvar.alphaSV*pi*rref.SV^2;
    % 
    %   else
    %       Rvar.RSV = (((Prmc.kr.SV*Lc.Leffsv)/((pi*rref.SV^2))^2)*(1-(Pvar.TMSV/(Prmc.kp.SV)))^(4/3))/nvar.nSV;
    %       Rvar.alphaSV = (1 - (Pvar.TMSV)/(Prmc.kp.SV))^(-2/3);
    %       Rvar.ASV = Rvar.alphaSV*pi*rref.SV^2;
    %   end
    % 
    % 
    % 
    % if Pvar.TMLV >= 0
    %     Rvar.RLV = (((Prmc.kr.LV*Lc.Lefflv)/((pi*rref.LV^2))^2)*(1+(Pvar.TMLV/(Prmc.kp.LV*Prmc.kl.LV)))^-4)/nvar.nLV;
    %     Rvar.alphaLV = ((Pvar.TMLV)/(Prmc.kp.LV*Prmc.kl.LV) + 1)^2;
    %     Rvar.ALV = Rvar.alphaLV*pi*rref.LV^2;
    % 
    % else
    %     Rvar.RLV = (((Prmc.kr.LV*Lc.Lefflv)/((pi*rref.LV^2))^2)*(1-(Pvar.TMLV/(Prmc.kp.LV)))^(4/3))/nvar.nLV;
    %     Rvar.alphaLV = (1 - (Pvar.TMLV)/(Prmc.kp.LV))^(-2/3);
    %     Rvar.ALV = Rvar.alphaLV*pi*rref.LV^2;
    % 
    % end


    % % %Poiseuille only
        Rvar.RSV = (((Prmc.kr.SV*Lc.Leffsv)/((pi*rref.SV^2))^2)*(1+(Pvar.TMSV/(Prmc.kp.SV*Prmc.kl.SV)))^-4)/nvar.nSV;
        Rvar.alphaSV = ((Pvar.TMSV)/(Prmc.kp.SV*Prmc.kl.SV) + 1)^2;
        Rvar.ASV = Rvar.alphaSV*pi*rref.SV^2;

        Rvar.RLV = (((Prmc.kr.LV*Lc.Lefflv)/((pi*rref.LV^2))^2)*(1+(Pvar.TMLV/(Prmc.kp.LV*Prmc.kl.LV)))^-4)/nvar.nLV;
        Rvar.alphaLV = ((Pvar.TMLV)/(Prmc.kp.LV*Prmc.kl.LV) + 1)^2;
        Rvar.ALV = Rvar.alphaLV*pi*rref.LV^2;


    % Total resistance
    Rvar.total = Rvar.RLA + Rvar.RSA + Rvar.RC + Rvar.RSV + Rvar.RLV;

    % Resistance convergence calculation
    mydiff = sqrt((Rvarprev.RSV-Rvar.RSV)^2 + (Rvarprev.RLV-Rvar.RLV)^2)/min(Rvar.RSV,Rvar.RLV); 

    % update venule wall shear stress calculations
    tauvar.tauSV = 32*muc.muSV*Qvar.QSV/(pi*Dvar.DSV^3);
    tauvar.tauLV = 32*muc.muLV*Qvar.QLV/(pi*Dvar.DLV^3);
  

    end
    
end
