function [z,timedata,spatialdata] = DAall(t,u,Pc,Dc,nc,Lc,muc,Prmc,Rc,Cc,Qc,...
  M0in,Lengthconstant,multiplier,respswitches,...
  krogh,formatvars,condrespsignalc,tauc, ind, rref,IOPin,Pvar)



%  Store new diameters in a structure
%  Initialize/store non-changing diameters
Dvar = Dc;
nvar = nc;
%  Store the changing diameters (again notice importance of consistent
%  field names
Dvar.DLA = u(1); zswitch(1) = 1;
Dvar.DSA = u(3); zswitch(3) = 1;

%  Go ahead and store activations for later use as well
Avar.ALA = u(2); zswitch(2) = 1;
Avar.ASA = u(4); zswitch(4) = 1;


%  Below has been moved to a separate m-file (as we do it a couple of times
%  in different locations)
%%% Essentially does the same thing as SetCState where it calculates flow,
%%% resistance, pressure, etc.

% if ind==
    [Pvar,Qvar,tauvar,Rvar,Cvar] = obtainnewflowquantselev(Dvar,Dc,...
        Pc,nvar,Lc,muc,Rc,Cc,Qc, Prmc,ind, rref,IOPin,Pvar);
    % [Pvar,Qvar,tauvar,Rvar,Cvar] = obtainnewflowquantscontrol(Dvar,Dc,...
        % Pc,nvar,Lc,muc,Rc,Cc,Qc, Prmc,ind, rref,IOPin);

% else 
    % [Pvar,Qvar,tauvar,Rvar,Cvar] = obtainnewflowquantscontrol(Dvar,Dc,Pc,nvar,Lc,muc,Rc,Cc,Qc, Prmc,ind, rref);
% end 

    timedata.flow.P = Pvar;
    timedata.flow.Q = Qvar;
    timedata.flow.tau = tauvar;
    timedata.flow.R = Rvar;
    timedata.flow.C = Cvar;
 

%%% Essentially the same as getscJAcontrol, but for after the control state
[condrespsignal,constotncu,spatialdata] = ...
  getscMRSRCR(formatvars,Dvar,Dc,Qvar,nvar,nc,Lc,M0in,...
  Lengthconstant,Prmc,krogh,multiplier);   


%  Store conducted response
timedata.meta.scr.LA = condrespsignal(1);         
timedata.meta.scr.SA = condrespsignal(2);

%  Store other metabolic info
%%% Create a vector for each junction in the vessel
jnctnames = {'ULA','LASA','SAC','CSV','SVLV',...
  'LVV'};
for fnc = 1:numel(jnctnames)
  fn = jnctnames{fnc};
  timedata.meta.sat.(fn) = spatialdata.endpts.O2Hbsat(fnc);
  timedata.meta.PO2.(fn) = spatialdata.endpts.PO2(fnc);
  timedata.meta.conc.(fn) = spatialdata.endpts.ATPconc(fnc);
  timedata.meta.CO2.(fn) = spatialdata.endpts.CO2(fnc);

end

vessnames = {'LA','SA','C','SV','LV'};
qO2consmean = mean(spatialdata.grid.qO2cons,1)';
for fnc = 1:numel(vessnames)
  fn = vessnames{fnc};
  timedata.meta.sat.(fn) = mean(spatialdata.grid.O2Hbsat(:,fnc));
  timedata.meta.PO2.(fn) = mean(spatialdata.grid.PO2(:,fnc));
  timedata.meta.conc.(fn) = mean(spatialdata.grid.ATPconc(:,fnc));
  timedata.meta.CO2.(fn) = mean(spatialdata.grid.CO2(:,fnc));
  timedata.meta.cons.(fn) = qO2consmean(fnc);
end 

 tissue_CO2_content = CO2converter(timedata.flow.Q.QLA*nvar.nLA, timedata.meta.CO2.LV);
 tissue_CO2_value = content_to_CO2(tissue_CO2_content);


  %  The response levels
  %  myogenic
  timedata.resp.myo.LA = (-Prmc.LA.Cmyo*...
    ((Pvar.PLA-IOPin)*Dvar.DLA/2*respswitches.myo.LA+...
    (Pc.Pla-IOPin)*Dc.DLA/2*(~respswitches.myo.LA)));
  timedata.resp.myo.SA = (-Prmc.SA.Cmyo*...
    ((Pvar.PSA-IOPin)*Dvar.DSA/2*respswitches.myo.SA+...
    (Pc.Psa-IOPin)*Dc.DSA/2*(~respswitches.myo.SA)));

  %  shear
  timedata.resp.shear.LA = (Prmc.LA.Cshear*...
    (tauvar.tauLA*respswitches.shear.LA+...
    tauc.tauLA*(~respswitches.shear.LA)));
  timedata.resp.shear.SA = (Prmc.SA.Cshear*...
    (tauvar.tauSA*respswitches.shear.SA+...
    tauc.tauSA*(~respswitches.shear.SA)));

  %  meta
  timedata.resp.meta.LA = Prmc.LA.Cmeta*...
    (condrespsignal(1)*respswitches.meta.LA+...
    condrespsignalc(1)*(~respswitches.meta.LA));
  timedata.resp.meta.SA = Prmc.SA.Cmeta*...
    (condrespsignal(2)*respswitches.meta.SA+...
    condrespsignalc(2)*(~respswitches.meta.SA));

  %  CO2
  timedata.resp.CO2.LA = Prmc.LA.Cco2*...
    (tissue_CO2_value*respswitches.CO2.LA+...
    tissue_CO2_value*(~respswitches.CO2.LA));
  timedata.resp.CO2.SA = Prmc.SA.Cco2*...
    (tissue_CO2_value*respswitches.CO2.SA+...
    tissue_CO2_value*(~respswitches.CO2.SA));


  %  Total
  timedata.resp.total.LA = timedata.resp.myo.LA+...
    timedata.resp.shear.LA+timedata.resp.meta.LA;
  timedata.resp.total.SA = timedata.resp.myo.SA+...
    timedata.resp.shear.SA+timedata.resp.meta.SA;


z = [
  %  Dinst for LA
  1/Prmc.td*2/Pc.Pla*((Pvar.PLA-IOPin)*Dvar.DLA/2-...
  (Prmc.LA.Cpass*exp(Prmc.LA.Cppass*...
  ((Dvar.DLA/Prmc.LA.D0)-1))+...
  (Prmc.LA.Cact*exp(-1*...
  (((Dvar.DLA/Prmc.LA.D0)-Prmc.LA.Cpact)/...
  Prmc.LA.Cppact)^2))*Avar.ALA));...

  %  Activation for CLA
  1/Prmc.ta*(1/(1+exp(...
  timedata.resp.myo.LA+timedata.resp.shear.LA+timedata.resp.meta.LA...
  +timedata.resp.CO2.LA-Prmc.LA.Cpptone))-Avar.ALA);...


  %  Dinst for SA
  1/Prmc.td*2/Pc.Psa*((Pvar.PSA-IOPin)*Dvar.DSA/2-...
  (Prmc.SA.Cpass*exp(Prmc.SA.Cppass*...
  ((Dvar.DSA/Prmc.SA.D0)-1))+...
  (Prmc.SA.Cact*exp(-1*...
  (((Dvar.DSA/Prmc.SA.D0)-Prmc.SA.Cpact)/...
  Prmc.SA.Cppact)^2))*Avar.ASA));...

  %  Activation for SA
  1/Prmc.ta*(1/(1+exp(...
  timedata.resp.myo.SA+timedata.resp.shear.SA+timedata.resp.meta.SA+...
  timedata.resp.CO2.SA - Prmc.SA.Cpptone))-Avar.ASA);...

  
];