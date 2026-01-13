function [ATPconcend,ATPconc,x] =  upconccalcanalytic(formatvars,Prmc,...
  q,Q,D,Lbegin,Lend,Sin,Cin,CO2in)

%%% Calculate ATP values

%discrete calculation 
x = linspace(Lbegin,Lend,formatvars.nodeslength);


gamma = Prmc.kd*pi*D/((1-Prmc.Hd)*Q);
k2 = pi/4*D^2*Prmc.Ht*Prmc.R0/((1-Prmc.Hd)*Q)-...
  pi/4*D^2*Prmc.Ht*Prmc.R0*Prmc.R1*q*x(1)/((1-Prmc.Hd)*Q^2*Prmc.C0*Prmc.Hd)-...
  pi/4*D^2*Prmc.Ht*Prmc.R0*Prmc.R1*Sin/((1-Prmc.Hd)*Q);
k3 = pi/4*D^2*Prmc.Ht*Prmc.R0*Prmc.R1*q/((1-Prmc.Hd)*Q^2*Prmc.C0*Prmc.Hd);
k4 = pi/4*D^2*Prmc.Ht*Prmc.R0/((1-Prmc.Hd)*Q);

if Sin > 0
  ATPconc = k2/gamma+k3/gamma*x-k3/(gamma^2)+...
    exp(gamma*(x(1)-x))*(Cin-k2/gamma-k3/gamma*x(1)+k3/(gamma^2));
else
  ATPconc = k4/gamma+exp(gamma*(x(1)-x))*(Cin-k4/gamma);
end


ATPconcend = ATPconc(end);

ATPconc = ATPconc;