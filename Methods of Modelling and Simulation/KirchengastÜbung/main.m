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
step  = 0.1;

for i = 1:4

    [theta_euler(i,:), r_euler(i,:)] = f_euler(phi0(i), Hn, N0, theta0, begin, step, stop);
    [theta_rk(i,:), r_rk(i,:)]       = f_rungeKutta(phi0(i), Hn, N0, theta0, begin, step, stop);

end
%%
line = 2;
hold on 
plot(r_euler(line,:),theta_euler(line,:))
plot(r_rk(line,:),theta_rk(line,:))
hold off