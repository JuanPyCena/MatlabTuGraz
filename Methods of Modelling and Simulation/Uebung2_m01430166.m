%% UE Methods of Modelling and Simulation - Raeuber Beute Modell
% Sonnleitner Felix, 01430166
% 08.11.2017

%% Lotka und Voltera

clc;
clear all;
close all;
delta_t = 0.0001;
simtime = 100;
H0 = 150;
L0 = 12;
fzuH = 0.3;
fabH = 0.025;

fabL = 0.2;
fzuL = 0.0015;

H = H0;
L = L0;
maxtime = simtime/delta_t
time = linspace(0,simtime,maxtime);
for t = 1:maxtime
    resultH{t} = H;
    resultL{t} = L;
    
    dH = fzuH * H - fabH * (H * L);  
    dL = fzuL * (L  * H) - fabL * L; 
    
    H = H + delta_t * dH;
    L = L + delta_t * dL;
end

resultH_m = cell2mat(resultH);
resultL_m = cell2mat(resultL);

hold on
title('Anzahl von Luchsen und Hasen nach Lotka und Voltera')
plot(time,resultH_m)
plot(time,resultL_m)
xlabel('Zeit in Jahren')
ylabel('Anzahl in Stück')
legend('Anzahl der Hasen','Anzahl der Luchse')
hold off

figure
plot(resultH_m,resultL_m,'k')
title('Phasendiagramm - Luchs- und Hasenpopulation')
xlabel('Hasenpopulation')
ylabel('Luchspopulation')
hold on
plot(133.33,12,'ok')
text(133.33,12-0.2,'Fixpunkt (133.3, 12)','fontsize',8)
hold off
%% STROGATZ 
clc;
clear all;
close all;

delta_t = 0.001;
simtime = 500;
H0 = 10;
L0 = 2;
r1 = 0.04;
r2 = r1;
Lc = 10;
Hc = 50;

H = H0;
L = L0;

maxtime = simtime/delta_t
time = linspace(0,simtime,maxtime);
for t = 1:maxtime
    resultH{t} = H;
    resultL{t} = L;
    dH = r1 * (1 - (L/Lc)) * H;
    dL = -r2 * (1- (H/Hc)) * L;
    
    H = H + delta_t * dH;
    L = L + delta_t * dL;
end

resultH_m = cell2mat(resultH);
resultL_m = cell2mat(resultL);

hold on
title('Anzahl von Luchsen und Hasen nach STROGATZ')
plot(time,resultH_m)
plot(time,resultL_m)
xlabel('Zeit in Jahren')
ylabel('Anzahl in Stück')
legend('Anzahl der Hasen','Anzahl der Luchse')
hold off

figure
plot(resultH_m,resultL_m,'k')
title('Phasendiagramm - Luchs- und Hasenpopulation')
xlabel('Hasenpopulation')
ylabel('Luchspopulation')
hold on
plot(Hc,Lc,'ok')
text(Hc,Lc-0.2,'Fixpunkt (133.3, 12)','fontsize',8)
hold off





