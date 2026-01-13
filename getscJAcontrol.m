function [condrespsignalc,krogh,Prmc,systemstatec] = getscJAcontrol(formatvars,Dc,Qc,nc,...
  Lc,M0c,Lengthconstant,Prmc)

%  See initializeformatvars

% 
% % large arteriole
% DA1nc = D.Dfeedc;
% 
% % small arteriole 
% DA2nc =  D.Dmidc;
% 
% DA3nc = D.Dcolc;

[krogh,qO2potentialcons,Prmc] = findkroghconstants(Dc,nc,Lc,M0c,Prmc,[]);


%consumption/length in each compartment based on control M0
qO2potentialconsvec = [qO2potentialcons.LA,qO2potentialcons.SA,...
  qO2potentialcons.C];

    Lstart = 0;
    
    %  Compartment vectors 
    %%% Make vectors of flow, diameter, length and whatever V is
    Qnc = [Qc.QLA Qc.QSA Qc.QC Qc.QSV Qc.QLV];
    Q = Qnc;
    Vnc = [Qc.QLA/(pi*(Dc.DLA/2)^2) Qc.QSA/(pi*(Dc.DSA/2)^2) ...
      Qc.QC/(pi*(Dc.DC/2)^2) ...
      Qc.QSV/(pi*(Dc.DSV/2)^2) Qc.QLV/(pi*(Dc.DLV/2)^2)];
    Dnc = [Dc.DLA Dc.DSA Dc.DC Dc.DSV Dc.DLV];
    Lnc = [Lstart Lc.Leffla Lc.Leffsa Lc.Leffc Lc.Leffsv Lc.Lefflv];

        
    %  Special for plot-making
%     JBplotvars = load('JBplotvars');
    %  For plotting only so far.
%     Lnc2 = [Lstart JBplotvars.L.LCIA JBplotvars.L.LEF ...
%       L.LFEED L.LMID L.LCAP L.LvMID L.LvFEED];
   
    %  Initializing vectors

    %  Saturation
    systemstatec.endpts.O2Hbsat = zeros(length(Qnc)+1,1);
    systemstatec.endpts.O2Hbsat(1,:) = Prmc.SA1;
    %  Concentration
    systemstatec.endpts.ATPconc = zeros(length(Qnc)+1,1);
    systemstatec.endpts.ATPconc(1,:) = Prmc.CA1;
    %  CO2
    systemstatec.endpts.CO2 = zeros(length(Qnc)+1,1);
    systemstatec.endpts.CO2(1,:) = Prmc.CO2init;
    %  Partial pressure of oxygen
    systemstatec.endpts.PO2 = zeros(length(Qnc)+1,1);
    systemstatec.endpts.PO2(1,:) = Prmc.PbA1;
    
%     %  Oxygen lost along each segment/compartment?  This was never used in
%     %  the past and, I believe, just stayed as a zero vector
%     systemstate.segs.efflux = zeros(length(Qnc),1);
%     %  Not currently used
%     consumptionnc = zeros(length(Qnc)+1,1);

%     %  I don't know what this is used for...distances are stored in x's
%     systemstate.segs.dists = zeros(length(Qnc),1);    
    %{
    Not currently used
    runoutu = zeros(1);
    runoutla = zeros(1);
    runoutsa = zeros(1);
    runoutc = zeros(1);
    %}
    %  Initializations.
    systemstatec.grid.O2Hbsat = zeros(formatvars.nodeslength,length(Q),1);
    systemstatec.grid.ATPconc = zeros(formatvars.nodeslength,length(Q),1);
    systemstatec.grid.CO2 = zeros(formatvars.nodeslength,length(Q),1);
    systemstatec.grid.x = zeros(formatvars.nodeslength,length(Q),1);    
    systemstatec.grid.qO2cons = zeros(formatvars.nodeslength-1,length(Q),1);
%     %  Not currently used
%     xwhole = zeros(nodeslength, length(Q) + 1, 1);
    systemstatec.grid.PO2 = zeros(formatvars.nodeslength,length(Q),1);
%     systemstate.grid.efflux = zeros(formatvars.nodeslength,length(Q),1);
%     %  Not currently used
%     area = zeros(nodeslength-1, length(Q) - 2 , 1);

    Ltotal = cumsum(Lnc);
        
%%%%%%%%%%%%%%%%%
% % %Non-control state, constant consumption in each segment, with uptake

%loops through LA, SA, C compartments
  for j = 1:(length(Qnc)-2)
    
    %%analytic calculation; 100% extraction is allowed
    %         disp(compartments{j});
    % ??? I switched to using Ltotal to find the start and end points of
    % the segments.  Not sure if this is right.  Required changing
    % uptakecalcanalytic and upconccalcanalytic
    [systemstatec.endpts.O2Hbsat(j+1),systemstatec.endpts.ATPconc(j+1),...
      systemstatec.endpts.CO2(j+1),systemstatec.grid.O2Hbsat(:,j),...
      systemstatec.grid.ATPconc(:,j),systemstatec.grid.CO2(:,j),...
      systemstatec.grid.qO2cons(:,j),systemstatec.grid.x(:,j)] =...
      uptakecalcanalytic(formatvars,Prmc,qO2potentialconsvec(1,j),Qnc(1,j),...
      Dnc(1,j),Ltotal(1,j),Ltotal(1,j+1),systemstatec.endpts.O2Hbsat(j),...
      systemstatec.endpts.ATPconc(j),systemstatec.endpts.CO2(j));
   end
  

  % JB-This doesn't give total oxygen consumption in the case when we run
  % out of oxygen.  We try using our qO2cons instead (an estimate).  Note
  % that the original estimate was messed up due to lack of a 1/6000 factor
  % though perhaps this was taken into account whenever interpretted
  % elsewhere, we will see.
  
%   constotncu = (pi*(Rtc.COL^2-Radc.COL^2)*nc.nCOL*Lc.LCOL+...
%     pi*(Rtc.FEED^2-Radc.FEED^2)*nc.nFEED*Lc.LFEED+...
%     pi*(Rtc.MID^2-Radc.MID^2)*nc.nMID*Lc.LMID+...
%     pi*(Rtc.CAP^2-Radc.CAP^2)*nc.nCAP*Lc.LCAP)*(M0c/6000)/totalvolume;

  constotncu = sum(sum(diff(systemstatec.grid.x(:,1:3)).*...
    bsxfun(@times,[nc.nLA,nc.nSA,nc.nC],...
    systemstatec.grid.qO2cons(:,1:3))))/krogh.totalvolume;
  

  for j = length(Qnc)-1:length(Qnc)
    systemstatec.endpts.O2Hbsat(j+1) = systemstatec.endpts.O2Hbsat(j);
    systemstatec.grid.O2Hbsat(:,j) = systemstatec.endpts.O2Hbsat(j);
    systemstatec.endpts.CO2(j+1) = systemstatec.endpts.CO2(j);
    systemstatec.grid.CO2(:,j) = systemstatec.endpts.CO2(j);
    systemstatec.endpts.PO2(j+1)= systemstatec.endpts.PO2(j);
%     systemstate.segs.efflux(j) = systemstate.segs.efflux(j-1);
    
    %% analytic calculation, 100% extraction
    % %         disp(compartments{j});

    [systemstatec.endpts.ATPconc(j+1),systemstatec.grid.ATPconc(:,j),...
      systemstatec.grid.x(:,j)] = ...
      upconccalcanalytic(formatvars,Prmc,qO2potentialconsvec(1,end),...
      Qnc(1,j),Dnc(1,j),Ltotal(1,j),Ltotal(1,j+1),...
      systemstatec.endpts.O2Hbsat(j),systemstatec.endpts.ATPconc(j),...
      systemstatec.endpts.CO2(j));
    
  end
  
  

    saturationtotal = [systemstatec.grid.O2Hbsat(:,1);...
      systemstatec.grid.O2Hbsat(2:end,2);...
      systemstatec.grid.O2Hbsat(2:end,3);...
      systemstatec.grid.O2Hbsat(2:end,4);...
      systemstatec.grid.O2Hbsat(2:end,5)];
    concentrationtotal = [systemstatec.grid.ATPconc(:,1);...
      systemstatec.grid.ATPconc(2:end,2);...
      systemstatec.grid.ATPconc(2:end,3);...
      systemstatec.grid.ATPconc(2:end,4);...
      systemstatec.grid.ATPconc(2:end,5)]; 
    CO2total = [systemstatec.grid.CO2(:,1);...
      systemstatec.grid.CO2(2:end,2);...
      systemstatec.grid.CO2(2:end,3);...
      systemstatec.grid.CO2(2:end,4);...
      systemstatec.grid.CO2(2:end,5)];
    Pbtotal = [systemstatec.grid.PO2(:,1);...
      systemstatec.grid.PO2(2:end,2);...
      systemstatec.grid.PO2(2:end,3);...
      systemstatec.grid.PO2(2:end,4);...
      systemstatec.grid.PO2(2:end,5)]; 

xs = [systemstatec.grid.x(:,1);systemstatec.grid.x(2:end,2);...
  systemstatec.grid.x(2:end,3);systemstatec.grid.x(2:end,4);...
  systemstatec.grid.x(2:end,5)];


    

%  Go ahead and store the "control state" saturation at the end of the
%  vessels
Prmc.satendcapc = systemstatec.grid.O2Hbsat(end,3);
    
% Plots of full oxygen saturation and ATP concentration curves!!  3/15/06
%{
figure(3)
subplot(3,1,1);
my_colors = 'rbkgmc';
for j = 1:size(systemstatec.grid.x,2)
  %%% Plot Oxyhemoglobin saturation on first subplot
  plot(systemstatec.grid.x(:,j),systemstatec.grid.O2Hbsat(:,j),...
    my_colors(j),'Linewidth',2);
  hold on
end
%  Make a "colororder" matrix for coloring the lines like how they were
%  colored in the past.  Each row of matrix corresponds to a color, red
%  ('r'), blue ('b'), black ('k'), green ('g'), magenta ('m').  (See
%  ColorSpec)
plot(Ltotal,systemstatec.endpts.O2Hbsat,'k.','Markersize',20);
subplot(3,1,2);
for j = 1:size(systemstatec.grid.x,2)
  plot(systemstatec.grid.x(:,j),systemstatec.grid.ATPconc(:,j),...
    my_colors(j),'Linewidth',2);
  hold on
end
%%% Plot ATP concentration on second subplot
plot(Ltotal,systemstatec.endpts.ATPconc,'k.','Markersize',20);
subplot(3,1,3);
for j = 1:size(systemstatec.grid.x,2)
  plot(systemstatec.grid.x(:,j),systemstatec.grid.CO2(:,j),...
    my_colors(j),'Linewidth',2);
  hold on
end
%%% Plot C02 on third subplot
plot(Ltotal,systemstatec.endpts.CO2,'k.','Markersize',20);

subplot(3,1,1);
xlabel('Distance (cm)');
ylabel('Oxyhemoglobin Saturation')
% subplot(3,1,1),plot(linspace(Ltotal(6), Ltotal(7),formatvars.nodeslength),systemstate.grid.O2Hbsat(:,6,1),'c','Linewidth',2)
subplot(3,1,2);
xlabel('Distance (cm)');
ylabel('ATP Concentration (\muM)')
% subplot(3,1,2),plot(linspace(Ltotal(6), Ltotal(7),formatvars.nodeslength),systemstate.grid.ATPconc(:,6,1),'c','Linewidth',2)
subplot(3,1,3);
xlabel('Distance (cm)');
ylabel('CO2')
%}

% Plot saturation vs position
figure(3)
my_colors = 'rbkgmc';
for j = 1:size(systemstatec.grid.x,2)
  %%% Plot Oxyhemoglobin saturation on first subplot
  plot(systemstatec.grid.x(:,j),systemstatec.grid.O2Hbsat(:,j),...
    my_colors(j),'Linewidth',2);
  hold on
end
plot(Ltotal,systemstatec.endpts.O2Hbsat,'k.','Markersize',20);
xlabel('Distance (cm)');
ylabel('Oxyhemoglobin Saturation')

% Calculating perfusion and consumption

%%% Calculate perfusion from krogh values
%%% Eqn 16 in publication 3
perfusionc(1) = Qc.QLA*nc.nLA/krogh.totalvolume*6000; 


    xmid = (Ltotal(2:end)+Ltotal(1:end-1))./2;
    
    %  We only want the signal in the LA (1) and SA (2)
    for j = 1:2
      %  Between segment/compartment midpoint and endpoint of lower venule
      inds = (xs > xmid(j)) & (xs < Ltotal(end-1));
      condrespsignalc(j) = trapz(xs(inds),concentrationtotal(inds).*...
        exp(-(xs(inds)-xmid(j))/Lengthconstant));
    end
% 
% figure (50)
% load("O2consumption.mat");
% O2Hbsatavg;
% my_colors = 'rbkgmc';
% for j = 1:size(systemstatec.grid.x,2)
%   plot(systemstatec.grid.x(:,j),systemstatec.grid.O2Hbsat(:,j),...
%     my_colors(j),'Linewidth',2, 'LineStyle', '--');
%   plot(systemstatec1.grid.x(:,j),systemstatec1.grid.O2Hbsat(:,j),...
%     my_colors(j),'Linewidth',2);
%   plot(systemstatec2.grid.x(:,j),systemstatec2.grid.O2Hbsat(:,j),...
%     my_colors(j),'Linewidth',2);
%   hold on
% end
% plot(systemstatec1.grid.x(:,3),O2Hbsatavg,...
%     '-.black', 'Linewidth',1);
% %plot(Ltotal,systemstatec.endpts.O2Hbsat,'k.','Markersize',20);
% %plot(Ltotal,systemstatec2.endpts.O2Hbsat,'k.','Markersize',20);
% %  Make a "colororder" matrix for coloring the lines like how they were
% %  colored in the past.  Each row of matrix corresponds to a color, red
% %  ('r'), blue ('b'), black ('k'), green ('g'), magenta ('m').  (See
% %  ColorSpec)
% plot(Ltotal,systemstatec.endpts.O2Hbsat,'k.','Markersize',20);
% plot(Ltotal,systemstatec1.endpts.O2Hbsat,'k.','Markersize',20);
% plot(Ltotal,systemstatec2.endpts.O2Hbsat,'k.','Markersize',20);
% xlabel('Distance (cm)');
% ylabel('Oxyhemoglobin Saturation')
% 
% 
