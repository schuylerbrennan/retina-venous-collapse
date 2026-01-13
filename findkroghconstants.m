function [krogh,qO2potentialcons,Prm] = findkroghconstants(D,n,L,M0,Prm,krogh)


  %%% Radius of each section of the vessel
  krogh.Rad.LA = D.DLA/2;
  krogh.Rad.SA = D.DSA/2;
  krogh.Rad.C = D.DC/2;
  krogh.Rad.SV = D.DSV/2;
  krogh.Rad.LV = D.DLV/2;  
  
  %%% Calculates the volume of the vessel using Eqn 15 from publication 3
  %%% N = (n_c*L_c)/(SUM(pi*r_ti^2*n_i*L_i)) where r_ti = r_i + d in first
  %%% 3 sections and r_ti = r_i in venuous sections. 
  krogh.vesselvolume = n.nLA*L.Leffla*pi*krogh.Rad.LA^2+...
    n.nSA*L.Leffsa*pi*krogh.Rad.SA^2+...
    n.nC*L.Leffc*pi*krogh.Rad.C^2+...
    n.nSV*L.Leffsv*pi*krogh.Rad.SV^2+...
    n.nLV*L.Leffsv*pi*krogh.Rad.LV^2;

  krogh.capvesselarea = n.nC*pi*krogh.Rad.C^2;

if ~isfield(krogh,'depth')
    
%   if ~isfield(Prm,'totalvolume') % first time around (haven't 
      %calculated total volume yet)

    krogh.aa = pi*(n.nLA*L.Leffla+n.nSA*L.Leffsa+ n.nC*L.Leffc);
    krogh.bb = 2*pi*(krogh.Rad.LA*n.nLA*L.Leffla+...
      krogh.Rad.SA*n.nSA*L.Leffsa+...
      krogh.Rad.C*n.nC*L.Leffc);
    krogh.cc = pi*(n.nLA*L.Leffla*(krogh.Rad.LA^2)+...
      n.nSA*L.Leffsa*(krogh.Rad.SA^2)+...
      n.nC*L.Leffc*(krogh.Rad.C^2)+...
      n.nSV*L.Leffsv*(krogh.Rad.SV^2)+...
      n.nLV*L.Lefflv*(krogh.Rad.LV^2))-...
      ((n.nC*L.Leffc)/(Prm.maxCAPdens*100));

    %%% Quadratic formula on krogh a b and c to calculate krogh depth (d)
    krogh.depth = (-krogh.bb+sqrt(krogh.bb^2-4*krogh.aa*krogh.cc))/...
      (2*krogh.aa);
%   else
      
    %  Two options for ctv = current_total_volume
    %  1.  Total volume of tissue+vessel conserved
%     ctv = Prm.totalvolume;
%     %  2.  Volume of tissue only conserved (vessel volume allowed to flux)
% %     ctv = krogh.vesselvolume+Prm.tissuevolume;
% 
%     Ln.LA = L.LLA*n.nLA;
%     Ln.SA = L.LSA*n.nSA;
%     Ln.C = L.LC*n.nC;
%     Ln.SV = L.LSV*n.nSV;
%     Ln.LV = L.LLV*n.nLV;
%     sum1 = Ln.C+Ln.LA+Ln.SA;
%     discriminant = -((((krogh.Rad.CAP-krogh.Rad.FEED)^2*Ln.CAP+...
%       (krogh.Rad.COL-krogh.Rad.FEED)^2*Ln.COL+...
%       (krogh.Rad.MID-krogh.Rad.FEED)^2*Ln.MID+...
%       krogh.Rad.vMID^2*Ln.vMID+krogh.Rad.vFEED^2*Ln.vFEED)*Ln.FEED+...
%       ((krogh.Rad.CAP-krogh.Rad.COL)^2*Ln.COL+...
%       (krogh.Rad.CAP-krogh.Rad.MID)^2*Ln.MID+...
%       krogh.Rad.vMID^2*Ln.vMID+krogh.Rad.vFEED^2*Ln.vFEED)*Ln.CAP+...
%       ((krogh.Rad.COL-krogh.Rad.MID)^2*Ln.MID+...
%       krogh.Rad.vMID^2*Ln.vMID+krogh.Rad.vFEED^2*Ln.vFEED)*Ln.COL+...
%       Ln.MID*(Ln.vFEED*krogh.Rad.vFEED^2+Ln.vMID*krogh.Rad.vMID^2))*pi-...
%       ctv*sum1)*pi;
% 
%     if discriminant < 0
%       error('Too many vessels...probably vessel volume > total volume');
%     end
%     krogh.depth = ((discriminant)^(1/2)+(-Ln.CAP*krogh.Rad.CAP-Ln.COL*...
%       krogh.Rad.COL-Ln.FEED*krogh.Rad.FEED-Ln.MID*krogh.Rad.MID)*pi)/...
%       (pi*sum1);
%   end
end

  %%% Calculate new radii for arterioles when accounting for krogh depth
  %%% (d). See eqn 15 in publication 3 for more detail.
  krogh.Rt.LA = krogh.Rad.LA+krogh.depth;
  krogh.Rt.SA = krogh.Rad.SA+krogh.depth;
  krogh.Rt.C = krogh.Rad.C+krogh.depth;
  krogh.Rt.SV = krogh.Rad.SV;
  krogh.Rt.LV = krogh.Rad.LV;

  % Rtissue = [Rtc.COL; Rtc.FEED; Rtc.MID; Rtc.CAP; Rtc.vMID; Rtc.vFEED];

  % totalvessellength = nc.nCOL*Lc.LCOL+nc.nFEED*Lc.LFEED+nc.nMID*Lc.LMID+...
  %   nc.nCAP*Lc.LCAP+nc.nvMID*Lc.LvMID+nc.nvFEED*Lc.LvFEED;
  %%% ACTUAL Eqn 15 from publication 3 that accounts for tissue region
  %%% width d (krogh depth)
  krogh.totalvolume = n.nLA*L.Leffla*pi*krogh.Rt.LA^2+...
    n.nSA*L.Leffsa*pi*krogh.Rt.SA^2+...
    n.nC*L.Leffc*pi*krogh.Rt.C^2+...
    n.nSV*L.Leffsv*pi*krogh.Rt.SV^2+...
    n.nLV*L.Lefflv*pi*krogh.Rt.LV^2;

  krogh.captotalarea = n.nC*pi*krogh.Rt.C^2;

  capdensityratio = krogh.capvesselarea/krogh.captotalarea;
  %%% Calculate how much volume tissue takes up given previous volume
  %%% calculations
  krogh.tissuevolume = krogh.totalvolume-krogh.vesselvolume;
  if ~isfield(Prm,'totalvolume')
    Prm.totalvolume = krogh.totalvolume;
    Prm.tissuevolume = krogh.tissuevolume;
  end
  krogh.volumefraction = krogh.vesselvolume/krogh.totalvolume;

  %%% Calculates capillary density by fully doing Eqn 15 from publication 3
  %%% instead of just denominator summation
  krogh.capdensity = n.nC*L.Leffc/krogh.totalvolume;
  

  %%% Oxygen consumption per vessel length for first 3 sections
  %%% Eqn 2 (q) from publication 3
  qO2potentialcons.LA = pi*M0/6000*(krogh.Rt.LA^2-krogh.Rad.LA^2);
  qO2potentialcons.SA = pi*M0/6000*(krogh.Rt.SA^2-krogh.Rad.SA^2);
  qO2potentialcons.C = pi*M0/6000*(krogh.Rt.C^2-krogh.Rad.C^2);

end