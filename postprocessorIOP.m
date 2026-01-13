function postprocessorIOP(varargin)

timedata = varargin{1};
spacedata = varargin{2};
IOP = varargin{3};

%D vs Pa LA and SA 
 for j = 1:length(IOP)
     DLA(j) = timedata{1,j}.var.D.LA(end);
     DSA(j) = timedata{1,j}.var.D.SA(end);
     ALA(j) = timedata{1,j}.var.A.LA(end);
     ASA(j) = timedata{1,j}.var.A.SA(end);
 end
 
 figure(27);
 hold on;
 plot(IOP/1333, DLA*10000, 'ro-')
 plot(IOP/1333, DSA*10000, 'bo-')
 legend('LA','SA')
 xlabel('IOP')
 ylabel('Diameter(micrometer)')
 title('Diameter vs IOP in LA and SA')
    box on

%  writematrix(DLA*10000, '1layerDLAPaIOP15.dat')
%  writematrix(DSA*10000, '1layerDSAPaIOP15.dat')

 
 %A vs Pa LA and SA
 figure(28);
 hold on;
 plot(IOP/1333, ALA, 'ro-')
 plot(IOP/1333, ASA, 'bo-')
 legend('LA','SA')
  xlabel('IOP')
 ylabel('Activation')
 title('Activation vs IOP in LA and SA')
   box on

