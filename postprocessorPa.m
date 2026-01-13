function postprocessorPa(varargin)
   timedata = varargin{1};
   spacedata = varargin{2};
   Pa = varargin{3};
   CpptoneLA = varargin{4};
   CpptoneSA = varargin{5};
   % elevated_IOP = varargin{6};
   muc = varargin{6};
   color = varargin{7};
   Prmc = varargin{8};
   IOPin = varargin{9};

%    close all;
   
%    mycolors = 'rgmcbyk'; 

n = numel(Pa);
count = 1;

% for i = 1:n
%     if (timedata1{1,i}.flow.P.Pcsv(end)-IOPin)/1333 > -.736
%         timedata{1,count} = timedata1{1,i};
%         spacedata{1,count} = spacedata1{1,i};
%         Pa(count) = Pa1(i);
% 
%         count = count + 1;
%     end
% end

   
   Pa = Pa./1333;
   

for colorcounter = 1:length(color)   

    % D vs t LA
  figure(20);
  hold on;
   for j = 1:length(Pa)
        %plot(timedata{1,j}.t, timedata{1, j}.var.D.LA*10000, mycolors(j),'Linewidth',1)
               plot(timedata{1,j}.t, timedata{1, j}.var.D.LA*10000, 'Linewidth',2)

   end
   legend(num2str(Pa(1,:)'));
   xlabel('time');
   ylabel('Diameter(micrometers)');
   title('Diameter vs Time in LA with various Pa values');
   box on


%D vs t SA
  figure(21);
  hold on;
   for j = 1:length(Pa)
        %plot(timedata{1,j}.t, timedata{1, j}.var.D.SA*10000, mycolors(j),'Linewidth',1)
               plot(timedata{1,j}.t, timedata{1, j}.var.D.SA*10000, 'Linewidth',2)

   end
   legend(num2str(Pa(1,:)'));
   xlabel('time');
   ylabel('Diameter(micrometers)');
   title('Diameter vs Time in SA with various Pa values');
   box on


%A vs t LA
 figure(22);
 hold on;
   for j = 1:length(Pa)
        %plot(timedata{1,j}.t, timedata{1, j}.var.A.LA, mycolors(j),'Linewidth',1)
              plot(timedata{1,j}.t, timedata{1, j}.var.A.LA, 'Linewidth',2)

   end
   legend(num2str(Pa(1,:)'));
   xlabel('time');
   ylabel('Activation');
   title('Activation vs Time in LA with various Pa values');
   axis([0 800 0 1]);
   box on

%A vs t SA
 figure(23);
 hold on;
   for j = 1:length(Pa)
        %plot(timedata{1,j}.t, timedata{1, j}.var.A.SA, mycolors(j),'Linewidth',1)
               plot(timedata{1,j}.t, timedata{1, j}.var.A.SA, 'Linewidth',2)

   end
   legend(num2str(Pa(1,:)'));
   xlabel('time');
   ylabel('Activation');
   title('Activation vs Time in SA with various M0 values');
   axis([0 800 0 1]);
   box on
   
%perfusion vs Pa 
% figure(24);
% 
% %  flow_tissue(iamhere,simc) = QSS.QTotalCalf/krogh.totalvolume*6000;
%  
% hold on;
%    for j = 1:length(Pa)
% %        perfusion(j) = timedata{1, j}.meta.perfusion.LA(end);
%        perfusion(j) = timedata{1, j}.flow.Q.QLA(end)*nc.nLA/krogh.totalvolume*6000;
%    end
%    plot(Pa, perfusion, 'o-');
%    xlabel('Pa')
%    ylabel('Perfusion')
%    title('Perfusion vs Pa')
%    axis([30 60 0 30]);
   
   
   
  
   
% Saturation vs x
figure(25);
hold on;
segmentcolors = 'bbbbbb';
   % segmentcolors = 'rbkgmc';
   for i = 1:length(Pa)
       j = 1;
       fn = fieldnames(spacedata{1,1}.grid.x);
       for fnc = 1:numel(fn)
           plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).O2Hbsat.(fn{fnc}),...
               segmentcolors(j),'Linewidth',1);
           j = j+1;
            
           
       end
       text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
           num2str(Pa(i)));
   end
    xlabel('Distance (cm)');
    ylabel('Saturation');
    title('Saturation across network for various Pa');
    legend((fn(:)'));
    axis([0 3 0 1]);
    box on

    
    
% PO2 vs x
 % figure(26);
 % hold on;
 %   segmentcolors = 'rbkgmc';
 %   for i = 1:length(Pa)
 %       j = 1;
 %       fn = fieldnames(spacedata{1,1}.grid.x);
 %       for fnc = 1:numel(fn)
 %           plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).PO2.(fn{fnc})/1333,...
 %               segmentcolors(j),'Linewidth',1);
 %           j = j+1;
 % 
 % 
 %       end
 %       text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).PO2.LV(end)/1333,...
 %           num2str(Pa(i)));
 %   end
 %    xlabel('Distance (cm)');
 %    ylabel('PO2 (mmHg)');
 %    title('PO2 across network for various Pa');
 %    legend((fn(:)'));
 %    axis([0 3 0 100]);
 %    box on
delTot = Pa*1333 - 15*1333;

%D vs Pa LA and SA 
 for j = 1:length(Pa)
     DLA(j) = timedata{1,j}.var.D.LA(end);
     DSA(j) = timedata{1,j}.var.D.SA(end);
     % DC(j) = ((32*muc.muC*timedata{1,j}.flow.Q.QC(end))/(pi*timedata{1,j}.flow.tau.tauC(end)))^(1/3);
     % DSV(j) = ((32*muc.muSV*timedata{1,j}.flow.Q.QSV(end))/(pi*timedata{1,j}.flow.tau.tauSV(end)))^(1/3);
     % DLV(j) = ((32*muc.muLV*timedata{1,j}.flow.Q.QLV(end))/(pi*timedata{1,j}.flow.tau.tauLV(end)))^(1/3);
     ALA(j) = timedata{1,j}.var.A.LA(end);
     ASA(j) = timedata{1,j}.var.A.SA(end);

     RLA(j) = timedata{1,j}.flow.R.RLA(end);
     RSA(j) = timedata{1,j}.flow.R.RSA(end);
     RC(j) = timedata{1,j}.flow.R.RC(end);
     RSV(j) = timedata{1,j}.flow.R.RSV(end);
     RLV(j) = timedata{1,j}.flow.R.RLV(end);

     AreaLV(j) = timedata{1,j}.flow.R.ALV(end);
     AreaSV(j) = timedata{1,j}.flow.R.ASV(end);
     alphaLV(j) = timedata{1,j}.flow.R.alphaLV(end);
     alphaSV(j) = timedata{1,j}.flow.R.alphaSV(end);
     alphaLA(j) = timedata{1,j}.flow.R.alphaLA(end);
     alphaSA(j) = timedata{1,j}.flow.R.alphaSA(end);

     f1(j) = timedata{1,j}.flow.P.f1(end);
     f2(j) = timedata{1,j}.flow.P.f2(end);
     f3(j) = timedata{1,j}.flow.P.f3(end);
     f4(j) = timedata{1,j}.flow.P.f4(end);
     f5(j) = timedata{1,j}.flow.P.f5(end);

     % E_vein(j) = timedata{1,j}.flow.R.E_vein(end);

     % delPla(j) = 2*(Pa(j)*1333-timedata{1,j}.flow.P.PLA(end));
     % delPsa(j) = 2*(Pa(j)*1333-delPla(j)-timedata{1,j}.flow.P.PSA(end));
     % f1(j) = delPla(j)/delTot(j);
     % f2(j) = delPsa(j)/delTot(j);
 end

 % 
 % if ~elevated_IOP
 %    PAc = Pa;
 %    ALAc = ALA;
 %    ASAc = ASA;
 %    DLAc = DLA;
 %    DSAc = DSA;
 %    save('controlAandD', 'PAc', 'ALAc', 'ASAc', 'DLAc', 'DSAc')
 % end
 % 
 % if elevated_IOP
 %    load controlAandD.mat;
 % end
 
 figure(27);
 hold on;
 if colorcounter == 1
 plot(Pa, DLA*10000, 'ro-')
 plot(Pa, DSA*10000, 'bo-')
 else
 plot(Pa, DLA*10000, 'mo')
 plot(Pa, DSA*10000, 'co-')
 end
     
 legend('LA','SA')
 xlabel('Pa')
 ylabel('Diameter(micrometer)')
 title('Diameter vs Pa in LA and SA')
    box on

%  writematrix(DLA*10000, '1layerDLAPaIOP15.dat')
%  writematrix(DSA*10000, '1layerDSAPaIOP15.dat')

 
 %A vs Pa LA and SA
 figure(28);
 hold on;
  if colorcounter == 1

 plot(Pa, ALA, 'ro-')
 plot(Pa, ASA, 'bo-')
  else
 plot(Pa, ALA, 'mo-')
 plot(Pa, ASA, 'co-')
  end

      
 
 legend('LA','SA')
  xlabel('Pa')
 ylabel('Activation')
 title('Activation vs Pa in LA and SA')
   box on




%Pa vs diameter for large arterioles
 figure(65);
 hold on;
 plot(Pa, DLA*10000, color(colorcounter))
 % plot(Pa, DLA*10000, 'ko')
  % plot(Pa, DLA*10000, 'bo-')
 %plot(Pa, DLAc*10000, 'ro-')
 legend('Control IOP','Elevated IOP')
 xlabel('Pa')
 ylabel('Diameter(micrometer)')
 title('Diameter vs Pa in LA')
    box on

%  writematrix(DLA*10000, '1layerDLAPaIOP15.dat')
%  writematrix(DSA*10000, '1layerDSAPaIOP15.dat')

%Pa vs diameter for small arterioles 
 figure(66);
 hold on;
 plot(Pa, DSA*10000, color(colorcounter))
 %plot(Pa, DSAc*10000, 'ro-')
 legend('Control IOP','Elevated IOP')
 xlabel('Pa')
 ylabel('Diameter(micrometer)')
 title('Diameter vs Pa in SA')
    box on


 
 %A vs Pa LA and SA
 figure(28);
 hold on;
 plot(Pa, ALA, 'ro-')
 plot(Pa, ASA, 'bo-')
 legend('LA','SA')
  xlabel('Pa')
 ylabel('Activation')
 title('Activation vs Pa in LA and SA')
   box on

%  writematrix(ALA, '1layerALAPaIOP15.dat')
%  writematrix(ASA, '1layerASAPaIOP15.dat')

 %A vs Pa in LA
 figure(67);
 hold on;
 plot(Pa, ALA, color(colorcounter))
 %plot(Pa, ALAc, 'ro-')
 legend('Control IOP','Elevated IOP')
  xlabel('Pa')
 ylabel('Activation')
 title('Activation vs Pa in LA')
   box on


 %A vs Pa in LA
 figure(68);
 hold on;
 plot(Pa, ASA, color(colorcounter))
 %plot(Pa, ASAc, 'ro-')
 legend('Control IOP','Elevated IOP')
  xlabel('Pa')
 ylabel('Activation')
 title('Activation vs Pa in SA')
   box on

 
% %P at nodes
% figure(29);
% hold on;
% L = [0, spacedata{1, 1}.grid.x.LA(end),spacedata{1, 1}.grid.x.SA(end),...
%     spacedata{1, 1}.grid.x.C(end), spacedata{1, 1}.grid.x.SV(end), spacedata{1, 1}.grid.x.LV(end)];
% fn = fieldnames(timedata{1, 1}.flow.P); 
% for j = 1:length(Pa)
%        P = zeros(length(Pa),numel(fn)+1);
%        P(j,1) = Pa(j);
%        for i = 1:numel(fn)
%            P(j,i+1) = timedata{1, j}.flow.P.(fn{i})(end)/1333;
%        end
% %       plot(L,P(j,:),mycolors(j),'Linewidth',1)
%              plot(L,P(j,:),'Linewidth',2)
% 
%              P;
% end
%     xlabel('Distance (cm)');
%     ylabel('Pressure');
%     title('Pressure at each node for all Pa');
%     legend(num2str(Pa(1,:)'));
%        box on

% SB 6/17 - saving activtion and diameter from getsJAcontrol to plot
% against CState6 


% %Resistance vs M0
%  figure(30);
%  hold on;
%     for i = 1:length(Pa)
%           LARes(i) = timedata{1, i}.flow.R.LA(end);
%           SARes(i) = timedata{1, i}.flow.R.SA(end);
%           CRes(i) = timedata{1, i}.flow.R.C(end);
%           SVRes(i) = timedata{1, i}.flow.R.SV(end);
%           LVRes(i) = timedata{1, i}.flow.R.LV(end);
%     end
% plot(Pa,LARes,'ro-');
% plot(Pa,SARes, 'bo-');
% plot(Pa,CRes,'ko-');
% plot(Pa,SVRes, 'go-');
% plot(Pa,LVRes,'mo-');
% title('Resistance vs Pa')
% xlabel('Pa')
% ylabel('Resistance')
% legend('LA', 'SA', 'C', 'SV', 'LV');

%Flow vs Pa 
figure(31);
hold on;
box on
    for i = 1:length(Pa)
          Flow(i) = timedata{1, i}.flow.Q.Qtotal(end);
    end
plot(Pa,Flow,color(colorcounter), 'LineWidth', 2);
 legend('Control IOP','Elevated IOP')
title('Total Flow vs Arterial Pressure')
xlabel('Arterial Pressure (mmHg)', 'FontSize', 18)
ylabel('Total Flow (cm^3/s)', 'FontSize', 18)
set(gca,'FontSize',18);
   box on

%writematrix(Flow, '1layerflowPaIOP25.dat')
% flowindex = find(Pa==40);
% Flow;
% 
% figure(32)
% hold on
% plot(Pa,Flow./Flow(flowindex),color(colorcounter), 'LineWidth', 2);
%  legend('Control IOP','Elevated IOP')
% title('Normalized Flow vs Arterial Pressure')
% xlabel('Arterial Pressure (mmHg)', 'FontSize', 18)
% ylabel('Normalized Flow', 'FontSize', 18)
% set(gca,'FontSize',18);
%    box on
   
%   keyboard
%     
% %Smeta at endpoints 
% figure(33);
% hold on;
% fn = fieldnames(timedata{1, 1}.meta.smeta); 
% for j = 1:length(Pa)
%        Smeta = zeros(length(Pa),numel(fn));
%        for i = 1:numel(fn)
%           Smeta(j,i) = timedata{1, j}.meta.smeta.(fn{i})(end);
%        end
% %       plot(L,Smeta(j,:),mycolors(j),'Linewidth',1)   
%             plot(L,Smeta(j,:),'Linewidth',2)    
% 
% end
%     xlabel('Distance (cm)');
%     ylabel('Smeta');
%     title('Smeta at each node for all Pa');
%     legend(num2str(Pa(1,:)'));
%     
%Components of Stone LA
 figure(34);
 hold on;
for j = 1:length(Pa)
    myoLA(j) = -timedata{1, j}.resp.myo.LA(end);
end
plot(Pa,myoLA,'r','Linewidth',1)

for j = 1:length(Pa)
    shearLA(j) = timedata{1, j}.resp.shear.LA(end);
end
plot(Pa,shearLA,'b','Linewidth',1)

for j = 1:length(Pa)
    metaLA(j) = timedata{1, j}.resp.meta.LA(end);
end
plot(Pa,metaLA,'g','Linewidth',1)

% for j = 1:length(Pa)
%     co2LA(j) = timedata{1, j}.resp.CO2.LA(end);
% end
% plot(Pa,co2LA,'k','Linewidth',1)
% plot(Pa,co2LA+metaLA,'c','Linewidth',1)
% plot(Pa, myoLA-shearLA-metaLA-co2LA+CpptoneLA,'m','Linewidth',1)
% plot(Pa, CpptoneLA.*ones(length(Pa),1),'k-')
  
xlabel('Pa')
ylabel('Component of S')
title('Large Arteriole Stone components')
legend('myo(+)', 'shear(-)', 'meta(-)') %, 'CO2(*)','meta+co2','Stone','C''tone')
box on

%Components of Stone SA
 figure(35);
 hold on;

for j = 1:length(Pa)
    myoSA(j) = -timedata{1, j}.resp.myo.SA(end);
end
plot(Pa,myoSA,'r','Linewidth',1)

for j = 1:length(Pa)
    shearSA(j) = timedata{1, j}.resp.shear.SA(end);
end
plot(Pa,shearSA,'b','Linewidth',1)

for j = 1:length(Pa)
    metaSA(j) = timedata{1, j}.resp.meta.SA(end);
end
plot(Pa,metaSA,'g','Linewidth',1)

% for j = 1:length(Pa)
%     co2SA(j) = timedata{1, j}.resp.CO2.SA(end);
% end
% plot(Pa,co2SA,'k','Linewidth',1)
% plot(Pa,co2SA+metaSA,'c','Linewidth',1)
% plot(Pa, myoSA-shearSA-metaSA-co2SA+CpptoneSA,'m','Linewidth',1)  
% plot(Pa, CpptoneSA.*ones(length(Pa),1),'k-')

xlabel('Pa')
ylabel('Component of S')
title('Small Arteriole Stone components')
legend('myo(+)', 'shear(-)', 'meta(-)') %,'co2(*)','meta+co2','Stone','C''tone')
box on
% 
% figure(36);
%  hold on;
% for j = 1:length(Pa)
%     myoLA(j) = -timedata{1, j}.resp.myo.LA(end);
% end
% plot(Pa,myoLA./myoLA(flowindex),'r','Linewidth',1)
% 
% for j = 1:length(Pa)
%     shearLA(j) = timedata{1, j}.resp.shear.LA(end);
% end
% plot(Pa,shearLA./shearLA(flowindex),'b','Linewidth',1)
% 
% for j = 1:length(Pa)
%     metaLA(j) = timedata{1, j}.resp.meta.LA(end);
% end
% plot(Pa,metaLA./metaLA(flowindex),'g','Linewidth',1)
% 
% % for j = 1:length(Pa)
% %     co2LA(j) = timedata{1, j}.resp.CO2.LA(end);
% % end
% % plot(Pa,co2LA./co2LA(flowindex),'k','Linewidth',1)
% % plot(Pa,co2LA+metaLA,'c','Linewidth',1)
% % plot(Pa, myoLA-shearLA-metaLA-co2LA+CpptoneLA,'m','Linewidth',2)
% 
% 
% xlabel('Pa')
% ylabel('Normmalized component of Stone')
% title('Large Arteriole Stone components')
% legend('myo(+)', 'shear(-)', 'meta(-)') %, 'CO2(*)')
% box on
% 
% %Components of Stone SA
%  figure(37);
%  hold on;
% 
% for j = 1:length(Pa)
%     myoSA(j) = -timedata{1, j}.resp.myo.SA(end);
% end
% plot(Pa,myoSA/myoSA(flowindex),'r','Linewidth',1)
% 
% for j = 1:length(Pa)
%     shearSA(j) = timedata{1, j}.resp.shear.SA(end);
% end
% plot(Pa,shearSA/shearSA(flowindex),'b','Linewidth',1)
% 
% for j = 1:length(Pa)
%     metaSA(j) = timedata{1, j}.resp.meta.SA(end);
% end
% plot(Pa,metaSA/metaSA(flowindex),'g','Linewidth',1)
% 
% % for j = 1:length(Pa)
% %     co2SA(j) = timedata{1, j}.resp.CO2.SA(end);
% % end
% % plot(Pa,co2SA/co2SA(flowindex),'k','Linewidth',1)
% % plot(Pa,co2SA+metaSA,'c','Linewidth',1)
% % plot(Pa, myoSA-shearSA-metaSA-co2SA+CpptoneSA,'m','Linewidth',2)  
% 
% xlabel('Pa')
% ylabel('Component of S')
% title('Small Arteriole Stone components')
% legend('myo(+)', 'shear(-)', 'meta(-)') %,'co2(*)')
% box on

% end
%     
%   
% 


%SB 6/17 - plotting diameter and activation solution from getsJAcontrol
%against diameter and activation from Cstate6 

 % 
 % 
 % figure(53);
 % hold on;
 % plot(Pa, DLA*10000, 'ro-')
 % plot(Pa, DSA*10000, 'bo-')
 % plot(PAc, DLAc*10000, 'r+-')
 % plot(PAc, DSAc*10000, 'b+-')
 % legend('LA elevated','SA elevated', 'LA control', 'SA control')
 % xlabel('Pa')
 % ylabel('Diameter(micrometer)')
 % title('Diameter vs Pa in LA and SA')
 %    box on

  

 % figure(52);
 % hold on;
 % plot(Pa, ALA, 'ro-')
 % plot(Pa, ASA, 'bo-')
 % plot(PAc, ALAc, 'r+-')
 % plot(PAc, ASAc, 'b+-')
 % yline(.5, '--black')
 % legend('LA elevated','SA elevated', 'LA control', 'SA control')
 % xlabel('Pa')
 % ylabel('Activation')
 % title('Activation vs Pa in LA and SA')
 %   box on


% figure(100)
%    for i = 1:length(Pa)
%        j = 1;
%        fn = fieldnames(spacedata{1,1}.grid.x);
%        for fnc = 1:numel(fn)
%            plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).O2Hbsat.(fn{fnc}),...
%                segmentcolors(j),'Linewidth',3);
%            j = j+1;
% 
%        end
%        text(spacedata{1, 1}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
%            num2str(Pa(i)),'FontSize', 16);
%         venous_satLV(i) = spacedata{1, i}.meta(2).O2Hbsat.LV(end);
%         venous_satSV(i) = spacedata{1, i}.meta(2).O2Hbsat.SV(end);
%         venous_satC(i) = spacedata{1, i}.meta(2).O2Hbsat.C(end);
%         venous_satSA(i) = spacedata{1, i}.meta(2).O2Hbsat.SA(end);
%         venous_satLA(i) = spacedata{1, i}.meta(2).O2Hbsat.LA(end);
%    end


%% Plotting resistances 

%Pa vs resistance for large arterioles
 figure(73);
 hold on;
 % plot(Pa, RLAc, 'b-', 'Linewidth', 2)
 xlabel('Pa')
 ylabel('Resistance')
 plot(Pa, RLA, color(colorcounter), 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs Pa in LA')
    box on

%  writematrix(DLA*10000, '1layerDLAPaIOP15.dat')
%  writematrix(DSA*10000, '1layerDSAPaIOP15.dat')

%Pa vs resistance for small arterioles 
 figure(74);
 hold on;
 % plot(Pa, RSAc, 'b-', 'Linewidth', 2)
 xlabel('Pa')
 ylabel('Resistance')
 plot(Pa, RSA, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')
  title('Resistance vs Pa in SA')
    box on


%Pa vs resistance for small venules
 figure(75);
 hold on;
 % plot(Pa, RSVc, 'b-', 'Linewidth', 2)
 xlabel('Pa')
 ylabel('Resistance')
 plot(Pa, RSV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs Pa in SV')
    box on

 %Pa vs resistance for large veins
 figure(76);
 hold on;
 % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
 xlabel('Pa')
 ylabel('Resistance')
 plot(Pa, RLV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs Pa in LV')
    box on


 %Pa vs resistance for capillaries
 figure(77);
 hold on;
 % plot(Pa, RCc, 'b-', 'Linewidth', 2)
 xlabel('Pa')
 ylabel('Resistance')
 plot(Pa, RC, color(colorcounter), 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs Pa in capillaries')
    box on




   figure (84)
hold on
   segmentcolors = 'cccgmygrkgmcrbkmcrgbkmcrgbkmcrgbk';
      % segmentcolors = 'bbbbbb';


   for i = 1:length(Pa)
       j = 1;
       nodes.Pa = Pa*1333;
       nodes.LASA = (timedata{1,i}.flow.P.PLA(end)-Pa*1333) + timedata{1,i}.flow.P.PLA(end);
       nodes.SAC = (timedata{1,i}.flow.P.PSA(end)-nodes.LASA) + timedata{1,i}.flow.P.PSA(end);
       nodes.CSV = (timedata{1,i}.flow.P.PC(end)-nodes.SAC) + timedata{1,i}.flow.P.PC(end);
       nodes.SVLV = (timedata{1,i}.flow.P.PSV(end)-nodes.CSV) + timedata{1,i}.flow.P.PSV(end);
       nodes.Pv = ones(1,length(Pa))*15*1333;

       % nodesc.Pa = Pa*1333;
       % nodesc.LASA = (timedatac{1,i}.flow.P.PLA(end)-Pa*1333) + timedatac{1,i}.flow.P.PLA(end);
       % nodesc.SAC = (timedatac{1,i}.flow.P.PSA(end)-nodesc.LASA) + timedatac{1,i}.flow.P.PSA(end);
       % nodesc.CSV = (timedatac{1,i}.flow.P.PC(end)-nodesc.SAC) + timedatac{1,i}.flow.P.PC(end);
       % nodesc.SVLV = (timedatac{1,i}.flow.P.PSV(end)-nodesc.CSV) + timedatac{1,i}.flow.P.PSV(end);
       nodesc.Pv = 1333*15;
%        nodes.LV = (timedata{1,i}.flow.P.PLV(end)-nodes.SA) + timedata{1,i}.flow.P.PSA(end);
       fn = fieldnames(spacedata{1,1}.grid.x); 
       pfn = fieldnames(nodes);
       

       for fnc = 1:numel(fn)
               yax = linspace(nodes.(pfn{fnc})(i), nodes.(pfn{fnc+1})(i), length(spacedata{1, i}.meta(2).x.(fn{fnc})));
               % yaxc = linspace(nodesc.(pfn{fnc}(i)), nodesc.(pfn{fnc+1})(i), length(spacedatac{1, i}.meta(2).x.(fn{fnc})));

%            test = ones(length(spacedata{1, i}.meta(2).x.(fn{fnc})))*timedata{1,i}.flow.P.(pfn{fnc+1})(end);
%            testc = ones(length(spacedatac{1, i}.meta(2).x.(fn{fnc})))*timedatac{1,i}.flow.P.(pfn{fnc+1})(end);
           plot(spacedata{1, i}.meta(2).x.(fn{fnc}),yax./1333,...
               segmentcolors(i),'Linewidth',2)
           yline(0)
           hold on
%            plot(spacedatac{1, i}.meta(2).x.(fn{fnc}),yaxc./1333,...
%                segmentcolors(j),'Linewidth',2)
           % j = j+1;
           clear yax


       end
       % text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
       %     num2str(Pa(i)));
   end
    xlabel('Distance (cm)');
    ylabel('Pressure');
    title('Pressure across network for various IOP');
%    axis([0 3 0 1]);
    box on

    %% Plotting pressure drop percentages (f1 and f2) vs Pa

    %Pa vs f1
    figure(85);
    hold on;
    xlabel('Pa')
    ylabel('f1')
    plot(Pa, f1, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f1 vs Pa')
    box on

    %Pa vs f2
    figure(86);
    hold on;
    xlabel('Pa')
    ylabel('f2')
    plot(Pa, f2, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f2 vs Pa')
    box on

    %Pa vs f3
    figure(87);
    hold on;
    xlabel('Pa')
    ylabel('f3')
    plot(Pa, f3, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f3 vs Pa')
    box on

    %Pa vs f4
    figure(88);
    hold on;
    xlabel('Pa')
    ylabel('f4')
    plot(Pa, f4, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f4 vs Pa')
    box on

    %Pa vs f5
    figure(89);
    hold on;
    xlabel('Pa')
    ylabel('f5')
    plot(Pa, f5, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f5 vs Pa')
    box on


    %% Plotting transmural pressure over system
    clear nodes
       figure (90)

hold on
   segmentcolors = 'cccgrbgrbmygrkgmcrbkmcrgbkmcrgbkmcrgbk';
      % segmentcolors = 'rrrrrr';


   for i = 1:length(Pa)
       j = 1;
       nodes.Pa = Pa(i)*1333 - IOPin;
       nodes.LASA = timedata{1,i}.flow.P.Plasa(end) - IOPin;
       nodes.SAC = timedata{1,i}.flow.P.Psac(end) - IOPin;
       nodes.CSV = timedata{1,i}.flow.P.Pcsv(end) - IOPin;
       nodes.SVLV = timedata{1,i}.flow.P.Psvlv(end) - IOPin;
       % nodes.Pv = ones(1,length(Pa))*15*1333 - IOPin;
       nodes.Pv = 15*1333 - IOPin;


       % nodesc.Pa = Pa*1333;
       % nodesc.LASA = (timedatac{1,i}.flow.P.PLA(end)-Pa*1333) + timedatac{1,i}.flow.P.PLA(end);
       % nodesc.SAC = (timedatac{1,i}.flow.P.PSA(end)-nodesc.LASA) + timedatac{1,i}.flow.P.PSA(end);
       % nodesc.CSV = (timedatac{1,i}.flow.P.PC(end)-nodesc.SAC) + timedatac{1,i}.flow.P.PC(end);
       % nodesc.SVLV = (timedatac{1,i}.flow.P.PSV(end)-nodesc.CSV) + timedatac{1,i}.flow.P.PSV(end);
       nodesc.Pv = 1333*15;
%        nodes.LV = (timedata{1,i}.flow.P.PLV(end)-nodes.SA) + timedata{1,i}.flow.P.PSA(end);
       fn = fieldnames(spacedata{1,1}.grid.x); 
       pfn = fieldnames(nodes);
       

       for fnc = 1:numel(fn)
               yax = linspace(nodes.(pfn{fnc}), nodes.(pfn{fnc+1}), length(spacedata{1, i}.meta(2).x.(fn{fnc})));
               % yaxc = linspace(nodesc.(pfn{fnc}(i)), nodesc.(pfn{fnc+1})(i), length(spacedatac{1, i}.meta(2).x.(fn{fnc})));

%            test = ones(length(spacedata{1, i}.meta(2).x.(fn{fnc})))*timedata{1,i}.flow.P.(pfn{fnc+1})(end);
%            testc = ones(length(spacedatac{1, i}.meta(2).x.(fn{fnc})))*timedatac{1,i}.flow.P.(pfn{fnc+1})(end);
           plot(spacedata{1, i}.meta(2).x.(fn{fnc}), yax./1333,...
               'Color', segmentcolors(i), 'Linewidth',2)
           hold on
%            plot(spacedatac{1, i}.meta(2).x.(fn{fnc}),yaxc./1333,...
%                segmentcolors(j),'Linewidth',2)
           % j = j+1;
           clear yax


       end
       % text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
       %     num2str(Pa(i)));
   end
    xlabel('Distance (cm)');
    ylabel('Transmural Pressure');
    title('Transmural pressure across network for various IOP');
%    axis([0 3 0 1]);
    box on

    %Pa vs area for large veins
    figure(91);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('Pa')
    ylabel('Area')
    plot(Pa, AreaLV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('Area vs Pa in LV')
    box on

    %Pa vs area for small veins
    figure(92);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('Pa')
    ylabel('Area')
    plot(Pa, AreaSV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('Area vs Pa in SV')
    box on


        %Pa vs area for large veins
    figure(93);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('Pa')
    ylabel('alpha')
    plot(Pa, alphaLV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('alpha vs Pa in LV')
    box on

    %Pa vs area for small veins
    figure(94);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('Pa')
    ylabel('alpha')
    plot(Pa, alphaSV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('alpha vs Pa in SV')
    box on

    %Pa vs P(alpha) for large veins
    figure(95);
    hold on;
    xlabel('Pa')
    ylabel('alpha')
    plot(Pa, alphaLA,color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP', 'Elevated IOP')
    title('alpha vs Pa in LA')
    box on

    %Pa vs P(alpha) for small veins
    figure(96);
    hold on;
    xlabel('Pa')
    ylabel('alpha')
    plot(Pa, alphaSA, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP', 'Elevated IOP')
    title('alpha vs Pa in SA')
    box on


    for i = 1:length(Pa)
       venous_sat(i) = spacedata{1,i}.meta(2).O2Hbsat.LV(end);
   end

    figure(18);
    hold on;
    plot(Pa, venous_sat, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    xlabel('Pa')
    ylabel('Venous Saturation')
    title('Pa vs Saturation', 'FontSize', 17);


 %  %E_vein vs Pa 
 % figure(91);
 % hold on;
 % % plot(Pa, RCc, 'b-', 'Linewidth', 2)
 % xlabel('Pa')
 % ylabel('E_vein')
 % plot(Pa, E_vein, color(colorcounter), 'Linewidth', 2)
 % legend('Control IOP','Elevated IOP')
 % 
 % title('E_vein vs Pa')
 %    box on



% figure(70);
% hold on;
% plot(Pa,venous_satLV,'mo-','Linewidth',2)
% % plot(Pa,venous_satSV,'g','Linewidth',2)
% % plot(Pa,venous_satC,'ko-','Linewidth',2)
% % plot(Pa,venous_satSA/spacedata{1, length(Pa)}.meta(2).O2Hbsat.SA(end),'bo-','Linewidth',2)
% % plot(Pa,venous_satLA/spacedata{1, length(Pa)}.meta(2).O2Hbsat.LA(end),'ro-','Linewidth',2)
% %legend('LV', 'SV', 'C', 'SA', 'LA')
% xlabel('Pa')
% ylabel('Saturation')
% title('Pa vs Saturation', 'FontSize', 17);





end




