% Übung Kirchengast
% Felix Sonnleitner, 01430166
clc;
clear all;
close all;

Re     = 6371;           %km
z      = 650;            %km
theta0 = 0;              %Grad
phi0   = [0 30 60 89.9]; %Grad
N0     = 315;            %N-units
Hn     = 7;              %km
r0     = Re;             %km

begin = r0;
stop  = Re + z;
step  = 1;

for i = 1:numel(phi0)
    [theta_euler(i,:), r_euler(i,:)] = f_euler(phi0(i), Hn, N0, theta0, begin, step, stop);
    [theta_rk(i,:), r_rk(i,:)]       = f_rungeKutta(phi0(i), Hn, N0, theta0, begin, step, stop);
end

%% Plot

% polar plot
figure
circle = linspace(0,2*pi,500);
polarplot(circle,Re.*ones(size(circle))-Re, 'o')
title(['h = ',num2str(step)])
ax = gca;
ax.ThetaLim = [-30 30];
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
hold on
for line = 1:numel(phi0)
    polarplot(theta_euler(line,:), r_euler(line,:).'-Re, '--')
    polarplot(theta_rk(line,:), r_rk(line,:).'-Re)
end
polarplot(circle,z.*ones(size(circle)))
rlim([r_euler(1)-300-Re  r_euler(end)-Re]);
hold off
legendstring = {'R_E = r_0', ...
                '\phi_{0,E} = 0','\phi_{0,RK} = 0',...
                '\phi_{0,E} = 30','\phi_{0,RK} = 30',...
                '\phi_{0,E} = 60','\phi_{0,RK} = 60',...
                '\phi_{0,E} = 89.9','\phi_{0,RK} = 89.9',...
                'Satellitenorbit'};
legend(legendstring)

