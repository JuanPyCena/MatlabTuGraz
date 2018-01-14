%% Uebung Methoden der Modellierung und Simulation E03 (Steiner)
% Modellanpassung
% Sonnleitner Felix, 01430166
% 29.11.2017

%% Beispiel 1: Modellanpassung einer Geraden
clear all
close all
clc
y_obs = [2.2; 2.6; 3.6; 3.7; 6.1; 7.6; 8.3; 8.8; 9.4];
x = [1:9]';
sigma = 0.5;

%Punktsch?tzung
A = [ones(9,1),x];     %Designmatrix
a = inv(A'*A) * A' * y_obs;
y_mod = a(1) + a(2) * x;

%Intervallsch?tzung
C = inv(A' * A);    %Kovarianzmatrix
var_d = C(1,1);
var_k = C(2,2);
covar = C(1,2);
korrelation = covar / sqrt(var_d * var_k);

%Berechnung des Restfehlers, Guete der Anpassung
dy = ((y_obs - y_mod) ./ sigma).^2;
chi_obs = sum(dy);

figure(1)
hold on
plot(x,y_obs,'o')
plot(x,y_mod)
title('Modellanpassung einer Geraden')
xlabel('x')
ylabel('y')
legend('Messdaten','Modellgerade','location','northwest')
hold off

%% Beispiel 2 - Modellanpassung eines Polynoms zweiten Grades
% y = ax^2 + bx + c      

clear all
close all
clc

y_obs = [7.3; 3.1; 11.7; 11.9; 15.5; 31; 33.9; 47; 61.5; 71.9];
x = [0:9]';
sigma = 4.3;

%Punktsch?tzung
A =[ones(10,1),x,x.^2];
a = inv(A' * A) * A' * y_obs;
y_mod = a(3)*x.^2 + a(2)*x + a(1);

%Intervallsch?tzung
C = inv(A' * A);
for i = 1:3
    for j = 1:3
        R(i,j) = C(i,j) / sqrt(C(i,i) * C(j,j));        
    end
end

%Berechnung des Restfehlers, G¸te der Anpassung
dy = (y_obs - y_mod) ./ sigma;
error = y_obs - y_mod;

neg_dy = zeros(10,1);
neg_dy(find(error < 0)) = error(find(error < 0));
pos_dy = zeros(10,1);
pos_dy(find(error > 0)) = error(find(error > 0));
chi_obs = sum(dy.^2);
chi_alpha = 7 + 2 * sqrt(2*7);

figure(2)
hold on
plot(x,y_obs,'o')
%plot(x,y_mod)
errorbar(x,y_mod,neg_dy,pos_dy)
title('Modellanpassung eins Polynom zweiten Grades')
xlabel('x')
ylabel('y')
legend('Messdaten','angepasste Funktion','location','northwest')
hold off


%% Beispiel 3: Anpassung einer Geraden
% mit Fehler

clear all
close all
clc

y_obs = [0.033333, 0.021053, 0.017544, 0.015385, 0.012821, 0.012195, 0.011111, 0.010870, 0.010526, 0.010256]';
x = [0.002083, 0.001010, 0.000680, 0.000508, 0.000265, 0.000210, 0.000099, 0.000072, 0.000046, 0.000015]';
sigma = [0.003333, 0.001330, 0.000923, 0.000710, 0.000493, 0.000446, 0.000370, 0.000354, 0.000332, 0.000316]';

A = [ones(10,1) ./ sigma,x ./sigma];
a = inv(A'*A)*A'*(y_obs ./sigma);
y_mod = a(2) * x + a(1);
dy = y_mod - y_obs;
chi_obs = sum(dy.^2);

P = polyfit(x,y_obs,1);
yfit = P(1)*x+P(2);

figure(3)
hold on
plot(x,y_obs,'o')
% plot(x,y_mod)
errorbar(x,y_mod,[0;0;0;0;dy(5);0;dy(7:10)],[dy(1:4);0;dy(6);0;0;0;0])
plot(x,yfit,'b-.');
title('Anpassung einer Geraden - mit Fehler')
xlabel('x')
ylabel('y')
legend('Messdaten','angepasste Gerade','Polyfit 1. Ordnung', 'location','northwest')
hold off

%% Beispiel 3: Anpassung einer Geraden
%ohne Fehler
% clear all
% close all
% clc

y_obs = [0.033333; 0.021053; 0.017544; 0.015385; 0.012195; 0.011111; 0.010870; 0.010526; 0.010256];
x = [0.002083; 0.001010; 0.00068; 0.000508; 0.000210; 0.000099; 0.000072; 0.000046; 0.000015];
sigma = ones(length(x),1);%[0.003333; 0.001330; 0.000923; 0.000710; 0.000446; 0.00037; 0.000354; 0.000332; 0.000316];

A = [ones(9,1) ./ sigma,x ./sigma];
a1 = inv(A'*A)*A'*(y_obs ./sigma);
y_mod = a1(2) * x + a1(1);

dy = y_mod - y_obs;
neg_dy = zeros(9,1);
neg_dy(find(dy < 0)) = dy(find(dy < 0));
pos_dy = zeros(9,1);
pos_dy(find(dy > 0)) = dy(find(dy > 0));
chi_obs = sum(dy.^2);

P = polyfit(x,y_obs,1);
yfit = P(1)*x+P(2);

figure(4)
hold on
plot(x,y_obs,'o')
% plot(x,y_mod)
errorbar(x,y_mod,pos_dy,neg_dy)
plot(x,yfit,'b-.');
title('Anpassung einer Geraden - ohne Fehler')
xlabel('x')
ylabel('y')
legend('Messdaten','angepasste Gerade','Polyfit 1. Ordnung','location','northwest')
hold off

















