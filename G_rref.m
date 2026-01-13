function [rref, etac] = G_rref(Dc, Pc, Prmc)
% Function used to find the reference radii of all of the different vessels using 2(rref + η) = D,
% solved for rref, and then using the eta equation from Cassani 2016, equation 2.26.

rref.LA = fzero(@(x) (2*(1-Prmc.nu^2)*Pc.TMLA * x^2)/(Prmc.E_artery*Prmc.h.LA)...
    + 2*x - Dc.DLA, Dc.DLA/2); 
rref.SA = fzero(@(x) (2*(1-Prmc.nu^2)*Pc.TMSA * x^2)/(Prmc.E_artery*Prmc.h.SA)...
    + 2*x - Dc.DSA, Dc.DSA/2); 
rref.SV = fzero(@(x) (2*(1-Prmc.nu^2)*Pc.TMSV * x^2)/(Prmc.E_vein*Prmc.h.SV)...
    + 2*x - Dc.DSV, Dc.DSV/2); 
rref.LV = fzero(@(x) (2*(1-Prmc.nu^2)*Pc.TMLV * x^2)/(Prmc.E_vein*Prmc.h.LV)...
    + 2*x - Dc.DLV, Dc.DLV/2); 

rref.C = fzero(@(x) (2*(1-Prmc.nu^2)*Pc.TMC * x^2)/(Prmc.E_vein*Prmc.h.C)...
    + 2*x - Dc.DC, Dc.DC/2);

%solving for eta in the control state
etac.LA = (Dc.DLA/2) - rref.LA;
etac.SA = (Dc.DSA/2) - rref.SA;
etac.SV = (Dc.DSV/2) - rref.SV;
etac.LV = (Dc.DLV/2) - rref.LV;

etac.C = (Dc.DC/2) - rref.C;

% pLA = [(2*(1-Prmc.nu^2*Pc.TMLA))/(Prmc.E_artery*Prmc.h.LA) 2 -Dc.DLA];
% rref.LA1 = roots(pLA);
% rref.LA = rref.LA1(2,1);
% 
% pSA = [(2*(1-Prmc.nu^2*Pc.TMSA))/(Prmc.E_artery*Prmc.h.SA) 2 -Dc.DSA];
% rref.SA1 = roots(pSA);
% rref.SA = rref.SA1(2,1);
% 
% pC = [(2*(1-Prmc.nu^2*Pc.TMC))/(Prmc.E_vein*Prmc.h.C) 2 -Dc.DC];
% rref.C1 = roots(pC);
% rref.C = rref.C1(2,1);
% 
% pSV = [(2*(1-Prmc.nu^2*Pc.TMSV))/(Prmc.E_vein*Prmc.h.SV) 2 -Dc.DSV];
% rref.SV1 = roots(pSV);
% rref.SV = rref.SV1(2,1);
% 
% 
% pLV = [(2*(1-Prmc.nu^2*Pc.TMLV))/(Prmc.E_vein*Prmc.h.LV) 2 -Dc.DLV];
% rref.LV1 = roots(pLV);
% rref.LV = rref.LV1(2,1);
% 
% etac.LA = (Dc.DLA/2) - rref.LA;
% etac.SA = (Dc.DSA/2) - rref.SA;
% etac.SV = (Dc.DSV/2) - rref.SV;
% etac.LV = (Dc.DLV/2) - rref.LV;
% 
% etac.C = (Dc.DC/2) - rref.C;


end