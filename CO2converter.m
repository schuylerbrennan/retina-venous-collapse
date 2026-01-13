% function [tissue_CO2, tissue_CO2_sat] = CO2converter(vessel_flow,venous_CO2, venous_CO2sat)
function [tissue_CO2_sat] = CO2converter(vessel_flow,venous_CO2sat)

load Yedata

Q=[.65,1,2]';
Qnorm = 0.68; %cm^3/sec
flow = Qnorm.*Q;

diff=data(:,2)-data(:,4);
ratio = diff./data(:,2);
coeffs = polyfit(flow,ratio,1);


% figure
% plot(flow, ratio, 'k.','Markersize',16)
% hold on
% q = linspace(0,flow(end));
% plot(q, coeffs(1).*q + coeffs(2),'k','Linewidth',2)
% xlabel('flow (cm^3/sec)')
% ylabel('ratio')


ratio_CO2 = coeffs(1)*vessel_flow + coeffs(2);

% tissue_CO2 = venous_CO2 - venous_CO2*ratio_CO2;

tissue_CO2_sat = venous_CO2sat - venous_CO2sat*ratio_CO2;

% diffc=diff./Q;
