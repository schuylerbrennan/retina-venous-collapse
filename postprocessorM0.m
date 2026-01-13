function postprocessorM0(varargin)
   timedata = varargin{1};
   spacedata = varargin{2};
   M0 = varargin{3};
   muc = varargin{4};
   color = varargin{5};
   Pa = varargin{6};
   IOPin = varargin{7};

n = numel(M0);
count = 1;

% for i = 1:n
%     if (timedata1{1,i}.flow.P.Pcsv(end)-IOPin)/1333 > -.736
%         timedata{1,count} = timedata1{1,i};
%         spacedata{1,count} = spacedata1{1,i};
%         M0(count) = M01(i);
% 
%         count = count + 1;
%     end
% end

   Pa = Pa/1333;
   
   
   %mycolors = 'rgmcbyk'; 

   for colorcounter = 1:length(color)

          % D vs t LA
  figure(20);
  hold on;
   for j = 1:length(M0)
        %plot(timedata{1,j}.t, timedata{1, j}.var.D.LA*10000, mycolors(j),'Linewidth',1)
               plot(timedata{1,j}.t, timedata{1, j}.var.D.LA*10000, 'Linewidth',2)

   end
   legend(num2str(M0(1,:)'));
   xlabel('time');
   ylabel('Diameter(micrometers)');
   title('Diameter vs Time in LA with various M0 values');
   box on


%D vs t SA
  figure(21);
  hold on;
   for j = 1:length(M0)
        %plot(timedata{1,j}.t, timedata{1, j}.var.D.SA*10000, mycolors(j),'Linewidth',1)
               plot(timedata{1,j}.t, timedata{1, j}.var.D.SA*10000, 'Linewidth',2)

   end
   legend(num2str(M0(1,:)'));
   xlabel('time');
   ylabel('Diameter(micrometers)');
   title('Diameter vs Time in SA with various M0 values');
   box on


%A vs t LA
 figure(22);
 hold on;
   for j = 1:length(M0)
        %plot(timedata{1,j}.t, timedata{1, j}.var.A.LA, mycolors(j),'Linewidth',1)
              plot(timedata{1,j}.t, timedata{1, j}.var.A.LA, 'Linewidth',2)

   end
   legend(num2str(M0(1,:)'));
   xlabel('time');
   ylabel('Activation');
   title('Activation vs Time in LA with various M0 values');
   axis([0 800 0 1]);
   box on

%A vs t SA
 figure(23);
 hold on;
   for j = 1:length(M0)
        %plot(timedata{1,j}.t, timedata{1, j}.var.A.SA, mycolors(j),'Linewidth',1)
               plot(timedata{1,j}.t, timedata{1, j}.var.A.SA, 'Linewidth',2)

   end
   legend(num2str(M0(1,:)'));
   xlabel('time');
   ylabel('Activation');
   title('Activation vs Time in SA with various M0 values');
   axis([0 800 0 1]);
   box on


% Saturation vs x
figure(25);
hold on;
   % segmentcolors = 'rbkgmc';
   segmentcolors = 'rrrrrr';
   for i = 1:length(M0)
       j = 1;
       fn = fieldnames(spacedata{1,1}.grid.x);
       for fnc = 1:numel(fn)
           plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).O2Hbsat.(fn{fnc}),...
               segmentcolors(j),'Linewidth',1);
           j = j+1;
            
           
       end
       text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
           num2str(M0(i)));
   end
    xlabel('Distance (cm)');
    ylabel('Saturation');
    title('Saturation across network for various M0');
    legend((fn(:)'));
    axis([0 3 0 1]);
    box on

   
% %D vs t LA
% figure(2);
% hold on;
%    for j = 1:length(M0)
%        plot(timedata{1,j}.t, timedata{1, j}.var.D.LA*10000, mycolors(j),'Linewidth',1)
%    end
%    legend(num2str(M0(1,:)'));
%    xlabel('time');
%    ylabel('Diameter(micrometers)');
%    title('Diameter vs Time in LA with various M0 values');
%    
% %D vs t SA
% figure(3);
% hold on;
%    for j = 1:length(M0)
%        plot(timedata{1,j}.t, timedata{1, j}.var.D.SA*10000, mycolors(j),'Linewidth',1)
%    end
%    legend(num2str(M0(1,:)'));
%    xlabel('time');
%    ylabel('Diameter(micrometers)');
%    title('Diameter vs Time in SA with various M0 values');
%   
%    
% %A vs t LA
%  figure(4);
%  hold on;
%    for j = 1:length(M0)
%        plot(timedata{1,j}.t, timedata{1, j}.var.A.LA, mycolors(j),'Linewidth',1)
%         
%    end
%    legend(num2str(M0(1,:)'));
%    xlabel('time');
%    ylabel('Activation');
%    title('Activation vs Time in LA with various M0 values');
%    axis([0 800 0 1]);
%    
% %A vs t SA
% figure(5);
% hold on;
%    for j = 1:length(M0)
%        plot(timedata{1,j}.t, timedata{1, j}.var.A.SA, mycolors(j),'Linewidth',1)
%         
%    end
%    legend(num2str(M0(1,:)'));
%    xlabel('time');
%    ylabel('Activation');
%    title('Activation vs Time in SA with various M0 values');
%    axis([0 800 0 1]);
   
% %perfusion vs M0 
% figure(6);
% hold on;
%    for j = 1:length(M0)
%        perfusion(j) = timedata{1, j}.meta.perfusion.LA(end);
%    end
%    plot(M0, perfusion, 'o-');
%    xlabel('M0')
%    ylabel('Perfusion')
%    title('Perfusion vs M0')
%   
%    
% %perfusion vs consumption
% figure(7);
% hold on;   
%    for j = 1:length(M0)
%        constot(j) = timedata{1, j}.meta.constot.LA(end);
%    end
%    plot(constot*6000, perfusion, 'o-');
%    xlabel('Consumption');
%    ylabel('Perfusion');
%    title('Perfusion vs Consumption');

% 
% %Saturation vs x
% figure(8);
% hold on;
% box on
%    segmentcolors = 'rbkgmckkkkk';
%    for i = 1:length(M0)
%        % j = 1;
%        fn = fieldnames(spacedata{1,1}.grid.x); 
%        for fnc = 1:numel(fn)
%            plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).O2Hbsat.(fn{fnc}),...
%                segmentcolors(fnc),'Linewidth',3);
% 
%        end
%        text(spacedata{1, 1}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
%            num2str(M0(i)),'FontSize', 16);
%         venous_satLV(i) = spacedata{1, i}.meta(2).O2Hbsat.LV(end);
%         venous_satSV(i) = spacedata{1, i}.meta(2).O2Hbsat.SV(end);
%         venous_satC(i) = spacedata{1, i}.meta(2).O2Hbsat.C(end);
%         venous_satSA(i) = spacedata{1, i}.meta(2).O2Hbsat.SA(end);
%         venous_satLA(i) = spacedata{1, i}.meta(2).O2Hbsat.LA(end);
%    end
% 
%     xlabel('Distance (cm)')
%     ylabel('Saturation')
%     title('Saturation across Compartments for Various Tissue Oxygen Demand Values','FontSize', 17);
%     %legend((fn(:)'),'FontSize',18);
%     set(gca,'FontSize',18);


%     figure(25);
% hold on;
%    segmentcolors = 'rbkgmc';
%    for i = 1:length(Pa)
%        j = 1;
%        fn = fieldnames(spacedata{1,1}.grid.x);
%        for fnc = 1:numel(fn)
%            plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).O2Hbsat.(fn{fnc}),...
%                segmentcolors(j),'Linewidth',1);
%            j = j+1;
% 
% 
%        end
%        text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
%            num2str(Pa(i)));
%    end
%     xlabel('Distance (cm)');
%     ylabel('Saturation');
%     title('Saturation across network for various Pa');
%     legend((fn(:)'));
%     axis([0 3 0 1]);
%     box on

%writematrix(venous_sat, '1layersatM0.dat')

% figure(18);
% hold on;
% plot(M0,venous_satLV,'m','Linewidth',2)
% plot(M0,venous_satSV,'g','Linewidth',2)
% plot(M0,venous_satC,'k','Linewidth',2)
% plot(M0,venous_satSA/spacedata{1, length(M0)}.meta(2).O2Hbsat.SA(end),'bo-','Linewidth',2)
% plot(M0,venous_satLA/spacedata{1, length(M0)}.meta(2).O2Hbsat.LA(end),'ro-','Linewidth',2)
% legend('LV', 'SV', 'C', 'SA', 'LA')
% xlabel('M_0')
% ylabel('Saturation')
% title('M0 vs Saturation', 'FontSize', 17);
% 


% %PO2 vs x
% figure(9);
% hold on;
%    for i = 1:length(M0)
%        j = 1;
%        fn = fieldnames(spacedata{1,1}.grid.x);
%        for fnc = 1:numel(fn)
%            plot(spacedata{1, i}.meta(2).x.(fn{fnc}),spacedata{1, i}.meta(2).PO2.(fn{fnc})/1333,...
%                segmentcolors(j),'Linewidth',1);
%            j = j+1;
% 
%        end
%        text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).PO2.LV(end)/1333,...
%            num2str(M0(i)));
% 
%    end
% 
%     xlabel('Distance (cm)');
%     ylabel('PO2 (mmHg)');
%      title('PO2 across network for various M0');
%     legend((fn(:)'));
%     axis([0 3 0 100]);
    
% D vs M0 LA and SA 
 for j = 1:length(M0)
     DLA(j) = timedata{1,j}.var.D.LA(end);
     DSA(j) = timedata{1,j}.var.D.SA(end);
     DC(j) = ((32*muc.muC*timedata{1,j}.flow.Q.QC(end))/(pi*timedata{1,j}.flow.tau.tauC(end)))^(1/3);
     DSV(j) = ((32*muc.muSV*timedata{1,j}.flow.Q.QSV(end))/(pi*timedata{1,j}.flow.tau.tauSV(end)))^(1/3);
     DLV(j) = ((32*muc.muLV*timedata{1,j}.flow.Q.QLV(end))/(pi*timedata{1,j}.flow.tau.tauLV(end)))^(1/3);
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
 end
 

 figure(27);
 hold on;
 if colorcounter == 1
 plot(M0, DLA*10000, 'ro-')
 plot(M0, DSA*10000, 'bo-')
 else
 plot(M0, DLA*10000, 'mo')
 plot(M0, DSA*10000, 'co-')
 end
     
 legend('LA','SA')
 xlabel('M_0')
 ylabel('Diameter(micrometer)')
 title('Diameter vs M_0 in LA and SA')
    box on

  %A vs Pa LA and SA
 figure(28);
 hold on;
  if colorcounter == 1

 plot(M0, ALA, 'ro-')
 plot(M0, ASA, 'bo-')
  else
 plot(M0, ALA, 'mo-')
 plot(M0, ASA, 'co-')
  end
 legend('LA','SA')
  xlabel('M_0')
 ylabel('Activation')
 title('Activation vs M_0 in LA and SA')
   box on

 
%M0 vs diameter for large arterioles
 figure(65);
 hold on;
 plot(M0, DLA*10000, color(colorcounter), 'LineWidth',2)
 % plot(Pa, DLA*10000, 'ko')
  % plot(Pa, DLA*10000, 'bo-')
 %plot(Pa, DLAc*10000, 'ro-')
 legend('Control IOP','Elevated IOP')
 xlabel('M_0')
 ylabel('Diameter(micrometer)')
 title('Diameter vs M_0 in LA')
    box on


 %M0 vs diameter for small arterioles 
 figure(66);
 hold on;
 plot(M0, DSA*10000, color(colorcounter), 'LineWidth',2)
 %plot(Pa, DSAc*10000, 'ro-')
 legend('Control IOP','Elevated IOP')
 xlabel('M_0')
 ylabel('Diameter(micrometer)')
 title('Diameter vs M_0 in SA')
    box on


 
 %A vs M0 LA and SA
 figure(28);
 hold on;
 plot(M0, ALA, 'ro-', 'LineWidth',2)
 plot(M0, ASA, 'bo-', 'LineWidth',2)
 legend('LA','SA')
  xlabel('M_0')
 ylabel('Activation')
 title('Activation vs M_0 in LA and SA')
   box on

%  writematrix(ALA, '1layerALAPaIOP15.dat')
%  writematrix(ASA, '1layerASAPaIOP15.dat')

 %A vs Pa in LA
 figure(67);
 hold on;
 plot(M0, ALA, color(colorcounter), 'LineWidth',2)
 %plot(Pa, ALAc, 'ro-')
 legend('Control IOP','Elevated IOP')
  xlabel('M_0')
 ylabel('Activation')
 title('Activation vs M_0 in LA')
   box on


 %A vs M0 in LA
 figure(68);
 hold on;
 plot(M0, ASA, color(colorcounter), 'LineWidth',2)
 %plot(Pa, ASAc, 'ro-')
 legend('Control IOP','Elevated IOP')
  xlabel('M0')
 ylabel('Activation')
 title('Activation vs M0 in SA')
   box on

  %Flow vs Pa 
figure(31);
hold on;
box on
    for i = 1:length(M0)
          Flow(i) = timedata{1, i}.flow.Q.Qtotal(end);
    end
plot(M0,Flow,color(colorcounter), 'LineWidth', 2);
 legend('Control IOP','Elevated IOP')
title('Total Flow vs Oxygen Demand')
xlabel('Oxygen Demand', 'FontSize', 18)
ylabel('Total Flow (cm^3/s)', 'FontSize', 18)
set(gca,'FontSize',18);
   box on

 %% Plotting resistances 

%M0 vs resistance for large arterioles
 figure(73);
 hold on;
 % plot(Pa, RLAc, 'b-', 'Linewidth', 2)
 xlabel('M_0')
 ylabel('Resistance')
 plot(M0, RLA, color(colorcounter), 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs M_0 in LA')
    box on

%  writematrix(DLA*10000, '1layerDLAPaIOP15.dat')
%  writematrix(DSA*10000, '1layerDSAPaIOP15.dat')

%M0 vs resistance for small arterioles 
 figure(74);
 hold on;
 % plot(Pa, RSAc, 'b-', 'Linewidth', 2)
 xlabel('M_0')
 ylabel('Resistance')
 plot(M0, RSA, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')
  title('Resistance vs M_0 in SA')
    box on


%M0 vs resistance for small venules
 figure(75);
 hold on;
 % plot(Pa, RSVc, 'b-', 'Linewidth', 2)
 xlabel('M_0')
 ylabel('Resistance')
 plot(M0, RSV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs M_0 in SV')
    box on

 %M0 vs resistance for large veins
 figure(76);
 hold on;
 % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
 xlabel('M_0')
 ylabel('Resistance')
 plot(M0, RLV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs M_0 in LV')
    box on


 %Pa vs resistance for capillaries
 figure(77);
 hold on;
 % plot(Pa, RCc, 'b-', 'Linewidth', 2)
 xlabel('M_0')
 ylabel('Resistance')
 plot(M0, RC, color(colorcounter), 'Linewidth', 2)
 legend('Control IOP','Elevated IOP')

 title('Resistance vs M_0 in capillaries')
    box on


figure (84)
hold on
   segmentcolors = 'gkgmcrbkmcrgbkmcrgbkmcrgbkmcrgbkm';
      % segmentcolors = 'rrrrrr';


   for i = 1:length(M0)
       j = 1;
       nodes.Pa = Pa*1333;
       nodes.LASA = (timedata{1,i}.flow.P.PLA(end)-Pa*1333) + timedata{1,i}.flow.P.PLA(end);
       nodes.SAC = (timedata{1,i}.flow.P.PSA(end)-nodes.LASA) + timedata{1,i}.flow.P.PSA(end);
       nodes.CSV = (timedata{1,i}.flow.P.PC(end)-nodes.SAC) + timedata{1,i}.flow.P.PC(end);
       nodes.SVLV = (timedata{1,i}.flow.P.PSV(end)-nodes.CSV) + timedata{1,i}.flow.P.PSV(end);
       nodes.Pv = 1333*15;

       % nodesc.Pa = Pa*1333;
       % nodesc.LASA = (timedatac{1,i}.flow.P.PLA(end)-Pa*1333) + timedatac{1,i}.flow.P.PLA(end);
       % nodesc.SAC = (timedatac{1,i}.flow.P.PSA(end)-nodesc.LASA) + timedatac{1,i}.flow.P.PSA(end);
       % nodesc.CSV = (timedatac{1,i}.flow.P.PC(end)-nodesc.SAC) + timedatac{1,i}.flow.P.PC(end);
       % nodesc.SVLV = (timedatac{1,i}.flow.P.PSV(end)-nodesc.CSV) + timedatac{1,i}.flow.P.PSV(end);
       %nodesc.Pv = 1333*15;
%        nodes.LV = (timedata{1,i}.flow.P.PLV(end)-nodes.SA) + timedata{1,i}.flow.P.PSA(end);
       fn = fieldnames(spacedata{1,1}.grid.x); 
       pfn = fieldnames(nodes);

       for fnc = 1:numel(fn)
               yax = linspace(nodes.(pfn{fnc}), nodes.(pfn{fnc+1}), length(spacedata{1, i}.meta(2).x.(fn{fnc})));
               % yaxc = linspace(nodesc.(pfn{fnc}(i)), nodesc.(pfn{fnc+1})(i), length(spacedatac{1, i}.meta(2).x.(fn{fnc})));

%            test = ones(length(spacedata{1, i}.meta(2).x.(fn{fnc})))*timedata{1,i}.flow.P.(pfn{fnc+1})(end);
%            testc = ones(length(spacedatac{1, i}.meta(2).x.(fn{fnc})))*timedatac{1,i}.flow.P.(pfn{fnc+1})(end);
           plot(spacedata{1, i}.meta(2).x.(fn{fnc}),yax./1333,...
               segmentcolors(i),'Linewidth',2)
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

    %% Plotting pressure drop percentages (f1 and f2) vs M0

    %M0 vs f1
    figure(85);
    hold on;
    xlabel('M_0')
    ylabel('f1')
    plot(M0, f1, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f1 vs M_0')
    box on

    %M0 vs f2
    figure(86);
    hold on;
    xlabel('M_0')
    ylabel('f2')
    plot(M0, f2, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f2 vs M_0')
    box on

    %M0 vs f3
    figure(87);
    hold on;
    xlabel('M_0')
    ylabel('f3')
    plot(M0, f3, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f3 vs M_0')
    box on

    %M0 vs f4
    figure(88);
    hold on;
    xlabel('M_0')
    ylabel('f4')
    plot(M0, f4, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f4 vs M_0')
    box on

    %M0 vs f5
    figure(89);
    hold on;
    xlabel('M_0')
    ylabel('f5')
    plot(M0, f5, color(colorcounter), 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('f5 vs M_0')
    box on

        clear nodes
       figure (90)

hold on
   segmentcolors = 'gbmygrkgmcrbkmcrgbkmcrgbkmcrgbk';
      % segmentcolors = 'rrrrrr';


   for i = 1:length(M0)
       j = 1;
       nodes.Pa = Pa*1333 - IOPin;
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
    ylabel('Transmural Pressure');
    title('Transmural pressure across network for various IOP');
%    axis([0 3 0 1]);
    box on



    %M0 vs area for large veins
    figure(91);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('M_0')
    ylabel('Area')
    plot(M0, AreaLV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('Area vs M_0 in LV')
    box on

    %M0 vs area for small veins
    figure(92);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('M_0')
    ylabel('Area')
    plot(M0, AreaSV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('Area vs M0 in SV')
    box on


        %Pa vs area for large veins
    figure(93);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('M_0')
    ylabel('alpha')
    plot(M0, alphaLV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('alpha vs M_0 in LV')
    box on

    %Pa vs area for small veins
    figure(94);
    hold on;
    % plot(Pa, RLVc, 'b-', 'Linewidth', 2)
    xlabel('M_0')
    ylabel('alpha')
    plot(M0, alphaSV, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP','Elevated IOP')

    title('alpha vs M_0 in SV')
    box on

     %Pa vs alpha for large arteries
    figure(95);
    hold on;
    xlabel('M_0')
    ylabel('alpha')
    plot(M0, alphaLA,color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP', 'Elevated IOP')
    title('alpha vs M_0 in LA')
    box on

    %Pa vs alpha for small arteries
    figure(96);
    hold on;
    xlabel('M_0')
    ylabel('alpha')
    plot(M0, alphaSA, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    legend('Control IOP', 'Elevated IOP')
    title('alpha vs M_0 in SA')
    box on


    %Components of Stone LA
 figure(34);
 hold on;
for j = 1:length(M0)
    myoLA(j) = -timedata{1, j}.resp.myo.LA(end);
end
plot(M0,myoLA,'r','Linewidth',1)

for j = 1:length(M0)
    shearLA(j) = timedata{1, j}.resp.shear.LA(end);
end
plot(M0,shearLA,'b','Linewidth',1)

for j = 1:length(M0)
    metaLA(j) = timedata{1, j}.resp.meta.LA(end);
end
plot(M0,metaLA,'g','Linewidth',1)

% for j = 1:length(Pa)
%     co2LA(j) = timedata{1, j}.resp.CO2.LA(end);
% end
% plot(Pa,co2LA,'k','Linewidth',1)
% plot(Pa,co2LA+metaLA,'c','Linewidth',1)
% plot(Pa, myoLA-shearLA-metaLA-co2LA+CpptoneLA,'m','Linewidth',1)
% plot(Pa, CpptoneLA.*ones(length(Pa),1),'k-')
  
xlabel('M_0')
ylabel('Component of S')
title('Large Arteriole Stone components')
legend('myo(+)', 'shear(-)', 'meta(-)') %, 'CO2(*)','meta+co2','Stone','C''tone')
box on

%Components of Stone SA
 figure(35);
 hold on;

for j = 1:length(M0)
    myoSA(j) = -timedata{1, j}.resp.myo.SA(end);
end
plot(M0,myoSA,'r','Linewidth',1)

for j = 1:length(M0)
    shearSA(j) = timedata{1, j}.resp.shear.SA(end);
end
plot(M0,shearSA,'b','Linewidth',1)

for j = 1:length(M0)
    metaSA(j) = timedata{1, j}.resp.meta.SA(end);
end
plot(M0,metaSA,'g','Linewidth',1)

% for j = 1:length(Pa)
%     co2SA(j) = timedata{1, j}.resp.CO2.SA(end);
% end
% plot(Pa,co2SA,'k','Linewidth',1)
% plot(Pa,co2SA+metaSA,'c','Linewidth',1)
% plot(Pa, myoSA-shearSA-metaSA-co2SA+CpptoneSA,'m','Linewidth',1)  
% plot(Pa, CpptoneSA.*ones(length(Pa),1),'k-')

xlabel('M_0')
ylabel('Component of S')
title('Small Arteriole Stone components')
legend('myo(+)', 'shear(-)', 'meta(-)') %,'co2(*)','meta+co2','Stone','C''tone')
box on

  for i = 1:length(M0)
       venous_sat(i) = spacedata{1,i}.meta(2).O2Hbsat.LV(end);
   end

    figure(18);
    hold on;
    plot(M0, venous_sat, color(colorcounter), 'LineStyle', '-', 'Linewidth', 2)
    xlabel('M_0')
    ylabel('Venous Saturation')
    title('M_0 vs Saturation', 'FontSize', 17);

% figure(10);
% hold on;
%   plot(M0, DLA*10000, 'ro-')
%   plot(M0, DSA*10000, 'bo-')
%   legend('LA','SA')
%   xlabel('M0')
%   ylabel('Diameter(micrometer)')
%   title('Diameter vs M0 in LA and SA')
% 
% %  writematrix(DLA*10000, '1layerDLAM0IOP15.dat')
% %  writematrix(DSA*10000, '1layerDSAM0IOP15.dat')
% 
% 
% %A vs M0 LA and SA
% figure(11);
% hold on;
%   plot(M0, ALA, 'ro-')
%   plot(M0, ASA, 'bo-')
%   legend('LA','SA')
%   xlabel('M0')
%   ylabel('Activation')
%   title('Activation vs M0 in LA and SA')
% 
% %  writematrix(ALA, '1layerALAM0IOP15.dat')
% %  writematrix(ASA, '1layerASAM0IOP15.dat')
% 
% 
% 
% % %P at nodes
% % figure(12);
% % hold on;
% %     L = [0, spacedata{1, 1}.grid.x.LA(end),spacedata{1, 1}.grid.x.SA(end),...
% %          spacedata{1, 1}.grid.x.C(end), spacedata{1, 1}.grid.x.SV(end), spacedata{1, 1}.grid.x.LV(end)];
% %     fn = fieldnames(timedata{1, 1}.flow.P); 
% %     for j = 1:length(M0)
% %         P = zeros(length(M0),numel(fn)+1);
% %         P(:,1) = 40;
% %         for i = 1:numel(fn)
% %             P(j,i+1) = timedata{1, j}.flow.P.(fn{i})(end)/1333;
% %         end
% %         plot(L,P(j,:),mycolors(j),'Linewidth',1)
% %          
% %     end
% %     xlabel('Distance (cm)');
% %     ylabel('Pressure');
% %     title('Pressure at each node for all M0');
% %     legend(num2str(M0(1,:)'));
% 
% % %Resistance vs M0
% % figure(13);
% % hold on;
% %     for i = 1:length(M0)
% %         LARes(i) = timedata{1, i}.flow.R.LA(end);
% %         SARes(i) = timedata{1, i}.flow.R.SA(end);
% %         CRes(i) = timedata{1, i}.flow.R.C(end);
% %         SVRes(i) = timedata{1, i}.flow.R.SV(end);
% %         LVRes(i) = timedata{1, i}.flow.R.LV(end);
% %     end
% % plot(M0,LARes,'ro-');
% % plot(M0,SARes, 'bo-');
% % plot(M0,CRes,'ko-');
% % plot(M0,SVRes, 'go-');
% % plot(M0,LVRes,'mo-');
% % title('Resistance vs M0')
% % xlabel('M0')
% % ylabel('Resistance')
% % legend('LA', 'SA', 'C', 'SV', 'LV');
% 
% %Smeta at endpoints 
% % figure(14);
% % mycolors = 'rgb';
% % hold on;
% % box on
% %     fn = fieldnames(timedata{1, 1}.meta.smeta); 
% %     Smeta = zeros(length(M0),numel(fn));
% %     for j = 1:length(M0)
% %         for i = 1:numel(fn)
% %            Smeta(j,i) = timedata{1, j}.meta.smeta.(fn{i})(end);
% %         end
% %     end
% %     plot(L,Smeta(1,:),'ro-','Linewidth',2)
% %     plot(L,Smeta(2,:),'go-','Linewidth',2)
% %     plot(L,Smeta(3,:),'bo-','Linewidth',2)
% %     xlabel('Distance (cm)', 'FontSize', 25);
% %     ylabel('S_m_e_t_a', 'FontSize', 25);
% %     title('Smeta at each node for all M0');
% %     legend({'M_0 = 0.5','M_0 = 1','M_0 = 2'},'FontSize', 16);
% %     set(gca,'FontSize',20);
% %      text(.28,.017,'LA','FontSize', 18);
% %      text(.8,.021,'SA','FontSize', 18);
% %      text(1.22,.0235,'C','FontSize', 18);
% %      text(1.6,.02,'SV','FontSize', 18);
% %      text(2.2,.01,'LV','FontSize', 18);
% 
% %end
% %    
% % %Components of Stone LA
% % figure(15);
% % hold on;
% % for j = 1:length(M0)
% %     myoLA(j) = -timedata{1, j}.resp.myo.LA(end);
% % end
% % plot(M0,myoLA,'r','Linewidth',1)
% % 
% % for j = 1:length(M0)
% %     shearLA(j) = timedata{1, j}.resp.shear.LA(end);
% % end
% % plot(M0,shearLA,'b','Linewidth',1)
% % 
% % for j = 1:length(M0)
% %     metaLA(j) = timedata{1, j}.resp.meta.LA(end);
% % end
% % plot(M0,metaLA,'g','Linewidth',1)
% %   
% % xlabel('M0')
% % ylabel('Component of S')
% % title('Large Arteriole Stone components')
% % legend('myo(+)', 'shear(-)', 'meta(-)')
% % 
% % 
% % %Components of Stone SA
% % figure(16);
% % hold on;
% %   for j = 1:length(M0)
% %      myoSA(j) = -timedata{1, j}.resp.myo.SA(end);
% %   end
% %   plot(M0,myoSA,'r','Linewidth',1)
% % 
% %   for j = 1:length(M0)
% %      shearSA(j) = timedata{1, j}.resp.shear.SA(end);
% %   end
% %   plot(M0,shearSA,'b','Linewidth',1)
% % 
% %   for j = 1:length(M0)
% %      metaSA(j) = timedata{1, j}.resp.meta.SA(end);
% %   end
% %   plot(M0,metaSA,'g','Linewidth',1)
% %   
% %     xlabel('M0')
% %     ylabel('Component of S')
% %     title('Small Arteriole Stone components')
% %     legend('myo(+)', 'shear(-)', 'meta(-)')
% 
% %flow
% figure(17);
% hold on
% box on
%     for i = 1:length(M0)
%           Flow(i) = timedata{1, i}.flow.Q.Qtotal(end);
%     end
% 
%     %writematrix(Flow, '1layerflowM0IOP25.dat')
%     Flow;
% plot(M0,Flow,'ko-', 'LineWidth', 2);
% title('Total Flow vs Tissue Oxygen Demand')
% xlabel('Oxygen Demand (M_0)')
% ylabel('Total Flow (cm^3/s)')
% set(gca,'FontSize',18);
% %axis([0 3.5 0 12e-4]);
% 
% figure(171);
% hold on
% box on
%     for i = 1:length(M0)
%           Flow(i) = timedata{1, i}.flow.Q.Qtotal(end);
%     end
% 
%     %writematrix(Flow, '1layerflowM0IOP25.dat')
% title('Total Flow vs Tissue Oxygen Demand')
% xlabel('Oxygen Demand (M_0)')
% ylabel('Total Flow (cm^3/s)')
% %ylim([.0002 .0012])
% plot(M0,Flow, color, 'LineWidth', 2);
% %legend('Control IOP', 'Elevated IOP')
% set(gca,'FontSize',18);
% hold on
% 
% 
% figure(60);
% hold on
% box on
%           %FlowLAScale = timedata{1,length(M0)}.flow.Q.QLA(end);
%           %FlowSAScale = timedata{1,length(M0)}.flow.Q.QSA(end);
%           %FlowCScale = timedata{1,length(M0)}.flow.Q.QC(end);
%     for i = 1:length(M0)
%           FlowLA(i) = timedata{1, i}.flow.Q.QLA(end);
%           FlowSA(i) = timedata{1, i}.flow.Q.QSA(end);
%           FlowC(i) = timedata{1, i}.flow.Q.QC(end);
%     end
% 
%     %writematrix(Flow, '1layerflowM0IOP25.dat')
%     %Flow;
% plot(M0,FlowLA,'ro-', 'LineWidth', 2);
% plot(M0,FlowSA,'bo-', 'LineWidth', 2);
% plot(M0,FlowC,'ko-', 'LineWidth', 2);
% legend('LA/LV', 'SA/SV', 'C');
% title('Flow vs Tissue Oxygen Demand by Compartment')
% xlabel('Oxygen Demand (M_0)')
% ylabel('Flow (cm^3/s)')
% set(gca,'FontSize',18);
% 
% 
% % SB 6/25 - saving activtion and diameter from getsJAcontrol to plot
% % against CState6 
% % if ~elevated_IOP
% %     M0c = M0;
% %     ALAc = ALA;
% %     ASAc = ASA;
% %     DLAc = DLA;
% %     DSAc = DSA;
% %     save('controlAandDandM0', 'ALAc', 'ASAc', 'DLAc', 'DSAc','M0c')
% % end
% % 
% % 
% % %plots for when IOP is elevated
% % if elevated_IOP
% %     load controlAandDandM0.mat;
% 
%  figure(54);
%  hold on;
%  plot(M0, DLA*10000, 'ro-')
%  plot(M0, DSA*10000, 'bo-')
%  % plot(M0c, DLAc*10000, 'r+-')
%  % plot(M0c, DSAc*10000, 'b+-')
%  legend('LA elevated','SA elevated', 'LA control', 'SA control')
%  xlabel('M0')
%  ylabel('Diameter(micrometer)')
%  title('Diameter vs Pa in LA and SA')
%     box on
% 
% 
% 
%  figure(55);
%  hold on;
%  plot(M0, ALA, 'ro-')
%  plot(M0, ASA, 'bo-')
%  %plot(M0c, ALAc, 'r+-')
%  %plot(M0c, ASAc, 'b+-')
%  yline(.5, '--black')
%  legend('LA elevated','SA elevated', 'LA control', 'SA control')
%  xlabel('M0')
%  ylabel('Activation')
%  title('Activation vs Pa in LA and SA')
%    box on
% 
% %M0 vs diameter for large arterioles
%  figure(65);
%  hold on;
%  %if ind
%     plot(M0, DLA*10000, 'bo-')
%     %plot(M0, DLA*10000, 'ro-')
%  %end
%  xlabel('M0')
%  ylabel('Diameter(micrometer)')
%  ylim 
%  legend('Elevated IOP','Control IOP')
% 
%  title('Diameter vs M0 in LA')
%     box on
% 
% %M0 vs diameter for small arterioles 
%  figure(66);
%  hold on;
%  plot(M0, DSA*10000, 'b+-')
%  %plot(M0, DSAc*10000, 'ro-')
%  legend('Elevated IOP','Control IOP')
%  xlabel('M0')
%  ylabel('Diameter(micrometer)')
%  title('Diameter vs M0 in SA')
%     box on
% 
% 
% 
% %  writematrix(ALA, '1layerALAPaIOP15.dat')
% %  writematrix(ASA, '1layerASAPaIOP15.dat')
% 
%  %A vs M0 in LA
%  figure(67);
%  hold on;
%  plot(M0, ALA, 'b+-')
%  %plot(M0, ALAc, 'ro-')
%  legend('Elevated IOP','Control IOP')
%   xlabel('M0')
%  ylabel('Activation')
%  title('Activation vs M0 in LA')
%    box on
% 
% 
%  %A vs M0 in LA
%  figure(68);
%  hold on;
%  plot(M0, ASA, 'b+-')
%  %plot(M0, ASAc, 'ro-')
%  legend('Elevated IOP','Control IOP')
%   xlabel('M0')
%  ylabel('Activation')
%  title('Activation vs M0 in SA')
%    box on
% 
% 
% % figure (84)
% % hold on
% %    segmentcolors = 'rkgmcrbkmcrgbk';
% %       % segmentcolors = 'rrrrrr';
% %       Pa = 40;
% % 
% %    for i = 1:length(M0)
% %        j = 1;
% %        nodes.Pa = Pa*1333;
% %        nodes.LASA = (timedata{1,i}.flow.P.PLA(end)-Pa*1333) + timedata{1,i}.flow.P.PLA(end);
% %        nodes.SAC = (timedata{1,i}.flow.P.PSA(end)-nodes.LASA) + timedata{1,i}.flow.P.PSA(end);
% %        nodes.CSV = (timedata{1,i}.flow.P.PC(end)-nodes.SAC) + timedata{1,i}.flow.P.PC(end);
% %        nodes.SVLV = (timedata{1,i}.flow.P.PSV(end)-nodes.CSV) + timedata{1,i}.flow.P.PSV(end);
% %        nodes.Pv = 1333*15;
% % 
% %        % nodesc.Pa = Pa*1333;
% %        % nodesc.LASA = (timedatac{1,i}.flow.P.PLA(end)-Pa*1333) + timedatac{1,i}.flow.P.PLA(end);
% %        % nodesc.SAC = (timedatac{1,i}.flow.P.PSA(end)-nodesc.LASA) + timedatac{1,i}.flow.P.PSA(end);
% %        % nodesc.CSV = (timedatac{1,i}.flow.P.PC(end)-nodesc.SAC) + timedatac{1,i}.flow.P.PC(end);
% %        % nodesc.SVLV = (timedatac{1,i}.flow.P.PSV(end)-nodesc.CSV) + timedatac{1,i}.flow.P.PSV(end);
% %        nodesc.Pv = 1333*15;
% % %        nodes.LV = (timedata{1,i}.flow.P.PLV(end)-nodes.SA) + timedata{1,i}.flow.P.PSA(end);
% %        fn = fieldnames(spacedata{1,1}.grid.x); 
% %        pfn = fieldnames(nodes);
% % 
% %        for fnc = 1:numel(fn)-1
% %                yax = linspace(nodes.(pfn{fnc})(i), nodes.(pfn{fnc+1})(i), length(spacedata{1, i}.meta(2).x.(fn{fnc})));
% %                % yaxc = linspace(nodesc.(pfn{fnc}(i)), nodesc.(pfn{fnc+1})(i), length(spacedatac{1, i}.meta(2).x.(fn{fnc})));
% % 
% % %            test = ones(length(spacedata{1, i}.meta(2).x.(fn{fnc})))*timedata{1,i}.flow.P.(pfn{fnc+1})(end);
% % %            testc = ones(length(spacedatac{1, i}.meta(2).x.(fn{fnc})))*timedatac{1,i}.flow.P.(pfn{fnc+1})(end);
% %            plot(spacedata{1, i}.meta(2).x.(fn{fnc}),yax./1333,...
% %                segmentcolors(i),'Linewidth',2)
% %            hold on
% % %            plot(spacedatac{1, i}.meta(2).x.(fn{fnc}),yaxc./1333,...
% % %                segmentcolors(j),'Linewidth',2)
% %            % j = j+1;
% %            clear yax
% % 
% % 
% %        end
% %        % text(spacedata{1, i}.meta(2).x.LV(end)+.05,spacedata{1, i}.meta(2).O2Hbsat.LV(end),...
% %        %     num2str(Pa(i)));
% %    end
% %     xlabel('Distance (cm)');
% %     ylabel('Pressure');
% %     title('Pressure across network for various IOP');
% % %    axis([0 3 0 1]);
% %     box on
   end 

end






