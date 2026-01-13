function [tauc,Pc,muc,M0c,Prmc] = SetCState(capdensity)
%   SetCState(MyoFEEDin,MyoMIDin,SFEEDin,SMIDin,Cmeta)
%% CONTROL STATE

%  Introduced structures...easier to pass back and forth between
%  programs and much more EFficient than save-load method of passing
%  parameters between programs.
%  To understand the change for one of the structures, all pressures got
%  stored into 1 structure P with individual fields Pa, P1, etc.  To pass
%  all pressures between programs, one needs only pass the structure P back
%  and forth and the fields (the individual pressures) will go along with
%  the structure.

%Equivalent to ReadInCons
%Rest pressures(SI: dyne*s/cm^5)
%%% Passive control state pressure. Incoming pressure is 40 and we subtract
%%% IOP which is 15. It is the pressure drop across the whole vessel.
Pc.Ppassive = 40*1333; 

%  JB-not really sure why an M0 vector is dEFined here...it seems like only
%  control values should be dEFined in a program called "SetCState"
% M0 = [1 8.276666 12];
%%% Control state O2 demand
M0c = 1;
% 
% NDS = length(M0);


%Wall shear stress
%%% Assumed values for wall shear stress
tauc.tauLA = 30;
tauc.tauSA = 30;
tauc.tauC = 15;
tauc.tauSV = 10;
tauc.tauLV = 10;


%Viscosities(SI: P); need to make diameter-dependent
%%% Assumed values for blood viscosity
muc.muLA = 0.022844;
muc.muSA = 0.020616;
muc.muC = 0.100116;
muc.muSV = 0.020915;
muc.muLV = 0.024440;


%%% Assumed pressure coming into the LA
Pc.Pa = 40*1333;
%%% Assumed IOP
Pc.IOP = 15*1333;
Pc.elevIOP = 25*1333;
Pc.IOPval = Pc.IOP;
%%% Assumed fraction of the pressure drop that occurs in the ...
%%% LA
Pc.f1 = 0.25;  %assumption
%%% SA
Pc.f2 = 0.4;   %assumption
%%% Calculated from Poiseuille's Law and symmetry
%%% Eqn 14 on publication 3
%%% Uses shear stress and blood viscosity
Pc.c1 = (tauc.tauSV/tauc.tauSA)^(4/3)*(muc.muSA/muc.muSV)^(1/3);
Pc.c2 = (tauc.tauLV/tauc.tauLA)^(4/3)*(muc.muLA/muc.muLV)^(1/3);
%%% Total fraction of pressure drop adds to 1. So it is 1 minus the rest of
%%% the fraction of the pressure drops
Pc.f3 =  1 - (1+Pc.c2)*Pc.f1 - (1+Pc.c1)*Pc.f2;  %from total drops with symmetry

%Pressure drops
%%% Total pressure drop over entire vessel
Pc.pressuredrop = Pc.Pa - Pc.IOP;
%%%%%Pc.Pv = 5*1333;
%%%%%Pc.pressuredrop = Pc.Pa- Pc.IOP;
%%% Calculated fraction of pressure drop for each section of vessel
Pc.delPla = Pc.f1*Pc.pressuredrop;
Pc.delPsa = Pc.f2*Pc.pressuredrop;
Pc.delPc =  Pc.f3*Pc.pressuredrop;
%%% Calculated pressure drop from symmetry of vessel
%%% Eqn 14 on publication 3
Pc.delPsv = Pc.delPsa*(tauc.tauSV/tauc.tauSA)^(4/3)*(muc.muSA/muc.muSV)^(1/3);  %from symmetry
Pc.delPlv = Pc.delPla*(tauc.tauLV/tauc.tauLA)^(4/3)*(muc.muLA/muc.muLV)^(1/3);  %from symmetry


%Pressures at endpoints and midpoints
%%% Pressure calculated at each node and the midpoint of each section
%%% Pressure at the end of the vessel (Pv)
Pc.Pv = Pc.IOP;
%%% Pressure between each section is the pressure at the proceeding section
%%% plus the pressure drop over that section
Pc.Psvlv = Pc.Pv + Pc.delPlv;
Pc.Pcsv = Pc.Psvlv + Pc.delPsv;
Pc.Psac = Pc.Pcsv + Pc.delPc;
Pc.Plasa = Pc.Psac + Pc.delPsa;
%%% Pressure at the midpoint of a section is the avg of the pressures at
%%% each end of that section
Pc.Pla = 0.5*(Pc.Pa + Pc.Plasa);
Pc.Psa = 0.5*(Pc.Plasa + Pc.Psac);
Pc.Pc = 0.5*(Pc.Psac + Pc.Pcsv);
Pc.Psv = 0.5*(Pc.Pcsv + Pc.Psvlv);
Pc.Plv = 0.5*(Pc.Psvlv + Pc.Pv);


%  OLD: Parameters  from Table 1 of Publication 2
%  These parameters come from "Theoretical model of blood flow
%  autoregulation: roles of myogenic, shear-dependent, and metabolic
%  responses" by Brian E. Carlson, Julia C. Arciero, and Timothy W. Secomb

%Parameters are updated when subtracting IOP from Passive pressure of 40;
%fit to Jeppesen data, as in 2013 IOVS Paper; new Cact value based on Fry's
%power law fit for Cact parameters
%%% C values for Stone equation that are found from previous model and data
Prmc.LA.Cpass = 361.48;
Prmc.LA.Cppass = 53.688;
Prmc.LA.Cact = 2114.2;
Prmc.LA.Cpact =0.93256;
Prmc.LA.Cppact = 0.11379;
Prmc.LA.Cmyo = 0.0091843;
Prmc.LA.Cdptone = 3.281;
Prmc.LA.Cshear = 0.0258;
Prmc.LA.D0 = 2*Prmc.LA.Cpass/Pc.Ppassive;
Prmc.LA.Cmyoco2tone = 1.0117;
Prmc.LA.Cdpco2tone = 3.907;
Prmc.LA.Cco2 = 0.0008;
   
Prmc.SA.Cpass = 197.01;
Prmc.SA.Cppass = 17.569;
Prmc.SA.Cact = 3089.6;
Prmc.SA.Cpact = 1.0249;
Prmc.SA.Cppact = 0.20373;
Prmc.SA.Cmyo = 0.025495;
Prmc.SA.Cdptone = 4.6238;
Prmc.SA.Cshear = 0.0258;
Prmc.SA.D0 = 2*Prmc.SA.Cpass/Pc.Ppassive;
Prmc.SA.Cmyoco2tone = 0.93172;
Prmc.SA.Cdpco2tone = 3.907;
Prmc.SA.Cco2 = 0.00013144;

%  We include other parameter definitions
%%% Time constants governing the rates of passive diameter and activation
%%% changes
Prmc.td = 1;
Prmc.ta = 60;

%  Some additional parameters that we add to this file because we should
%  avoid multiple dEFinitions in different files (here these were dEFined
%  in both getscJAconrtol and getscMRSRCR)
%  parameters H1564 Table 1: Theoretical model of metabolic blood
% flow regulation:  roles of ATP release by red blood cells and conducted
% responses
%%% A further explanation of these parameters is shown in publication 3
%  Oxygen CAPacity of rbcs
Prmc.C0 = 0.5;
%  Discharge hematocrit
Prmc.Hd = 0.4;
%  Tube hematocrit
Prmc.Ht = 0.3;
%  Half-maximal Hb-hemoglobin saturation
Prmc.P50 = 26*1333;
%  Initial blood partial pressure of oxygen at entry to system...assume 100
%  mmHg
Prmc.PbA1 = 100*1333;
%  The power in the hill equation for saturation coEFficient (n_H)
Prmc.power = 2.7;
%initial saturation
%%% Use the Hill equation to calculate saturation
%%% S = P^n_b / (P^n_S0 + P^n_b)
%%% Prmc.power is n and Prmc.PbA1/Prmc.P50 is P_b
Prmc.SA1 = (Prmc.PbA1/Prmc.P50)^Prmc.power/(1+(Prmc.PbA1/Prmc.P50)^Prmc.power);
%  We dEFine the initial ATP concentration here
Prmc.CA1 = 0.5;
%  We dEFine the initial CO2 here
Prmc.CO2init = 0.5;

%  ATP release: release as a function of oxygen saturation = R(S(x)) =
%  R_0*(1-R_1*S(x))
Prmc.R0 = 1.4;
Prmc.R1 = 0.891;


%  ATP decay rate
Prmc.kd = 2e-4;

%  Maximum CAPillary density, used for krogh calculations?
Prmc.maxCAPdens = capdensity;

%  Provides a threshold below which we no longer transfer between
%  compartments.
Prmc.frac = 0.1;
%  By making this 0, the threshold basically gets eliminated.
% Prmc.frac = 1e-1;

%  Use a const (true) or variable (false) D0
Prmc.constD0 = true;

%  Maximum time to integrate to (in seconds) and inner loop NN (DAsolver)
Prmc.tmax = 1200;
Prmc.NNmax = 10;
Prmc.fulltime = false;

%  Control state activation levels
%%% Activation is spread evenly in LA and SA and only occurs there. Total
%%% of 1, so each section gets 0.5 Activation
Prmc.LA.ContAct = 0.5;
Prmc.SA.ContAct = 0.5;



end