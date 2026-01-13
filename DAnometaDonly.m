function response = DAnometaDonly(time,diameters,P,Prm)

%  diameters- LA, SA
%  P-pressure structure
%  Prm-parameter structure (for ode constants/coefficients)
%%% ODE for Diameter only, not Activation.
%%% Sovled for LA and SA
%%% dD/dt = (1/tau_d)(D_c/T_c)(T_i - T_total)

%%% 1/Prm.td = i/tau_d
%%% 2/P.Pla = D_c/T_c
%%% ((P.Pla-P.IOP)*diameters(1)/2 = T_i, from T = delP*D/2
%%% The rest is T_Total calculated from T_pass + AT^max_act with A = 0.5

response = [1/Prm.td*2/P.Pla*((P.Pla-P.IOP)*diameters(1)/2-...
  (Prm.LA.Cpass*exp(Prm.LA.Cppass*((diameters(1)/Prm.LA.D0)-1))+...
  (Prm.LA.Cact*exp(-1*(((diameters(1)/Prm.LA.D0)-Prm.LA.Cpact)/...
  Prm.LA.Cppact)^2))*Prm.LA.ContAct));...
  1/Prm.td*2/P.Psa*((P.Psa-P.IOP)*diameters(2)/2-...
  (Prm.SA.Cpass*exp(Prm.SA.Cppass*((diameters(2)/Prm.SA.D0)-1))+...
  (Prm.SA.Cact*exp(-1*(((diameters(2)/Prm.SA.D0)-Prm.SA.Cpact)/...
  Prm.SA.Cppact)^2))*Prm.SA.ContAct))];



