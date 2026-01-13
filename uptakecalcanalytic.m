%This function reads in values for C0, Prmc.Hd, Prmc.Ht, Prmc.R0, Prmc.R1, consumption, flow
%rate, mean flow velocity, beginning and end position points for particular
%vessel segment, and incoming saturation level in order to calculate
%saturation levels and ATP levels while stepping down the vessel (each
%segment is cut into nodeslength pieces using linspace).  This function
%assumes constant consumption along each vessel.  Function can be used for
%control and non-control states. This function assumes there is ATP uptake
%(given by ku) by the endothelial wall throughout pathway.

function [O2Hbsatend,ATPconcend,CO2end,O2Hbsat,ATPconc,CO2,qO2cons,x] = ...
  uptakecalcanalytic(formatvars,Prmc,q,Q,D,Lbegin,Lend,Sin,Cin,CO2in)

%%% Solve Eqns 2 and 5 from publication 3
%%% This relates to calculateing ATP concentration and Saturation

%discrete calculation 
x = linspace(Lbegin,Lend,formatvars.nodeslength)';


gamma = Prmc.kd*pi*D/((1-Prmc.Hd)*Q);
k2 = pi/4*D^2*Prmc.Ht*Prmc.R0/((1-Prmc.Hd)*Q)-...
  pi/4*D^2*Prmc.Ht*Prmc.R0*Prmc.R1*q*x(1)/...
  ((1-Prmc.Hd)*Q^2*Prmc.C0*Prmc.Hd)-...
  pi/4*D^2*Prmc.Ht*Prmc.R0*Prmc.R1*Sin/((1-Prmc.Hd)*Q);
k3 = pi/4*D^2*Prmc.Ht*Prmc.R0*Prmc.R1*q/((1-Prmc.Hd)*Q^2*Prmc.C0*Prmc.Hd);
k4 = pi/4*D^2*Prmc.Ht*Prmc.R0/((1-Prmc.Hd)*Q);
k5 = pi/4*D^2*Prmc.Ht*Prmc.R0*(1-0.1*Prmc.R1)/((1-Prmc.Hd)*Q);


qO2cons = q*ones(numel(x)-1,1);
O2Hbsat = q*(x(1)-x)./(Q*Prmc.C0*Prmc.Hd)+Sin;
CO2 = -0.81*q*(x(1)-x)/(Q*Prmc.C0*Prmc.Hd)+CO2in;
ATPconc = k2/gamma+k3/gamma*x-k3/(gamma^2)+...
  exp(gamma*(x(1)-x))*(Cin-k2/gamma-k3/gamma*x(1)+k3/(gamma^2));
inds = (O2Hbsat <= 0);
if any(inds)
  O2Hbsat(inds) = 0;
  qO2cons(inds(1:end-1)) = 0;
  tmp = find(inds); O2Hbsat0ind = tmp(1);
  %CO2(O2Hbsat0ind:end) = -0.81*qO2cons(end)*(x(1)-x(O2Hbsat0ind:end))/(Q*Prmc.C0*Prmc.Hd)+CO2in;
  ATPconc(O2Hbsat0ind:end) = k4/gamma+exp(gamma*(x(O2Hbsat0ind)-...
    x(O2Hbsat0ind:end)))*(ATPconc(O2Hbsat0ind)-k4/gamma);
end



O2Hbsatend = O2Hbsat(end);
ATPconcend = ATPconc(end);
CO2end = CO2(end);