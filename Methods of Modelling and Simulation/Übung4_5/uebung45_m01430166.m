% Uebung E04/E05 Method of Modelling and Simulation
% Sonnleitner Felix, 01430166
% 09.12.2017

clear all
close all
clc

%% B1: Anpassung von Polynomen verschiedener Ordnung
t_obs = [-0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9]';
y_obs = [81 50 35 27 26 60 106 189 318 520]';
sigma = sqrt(y_obs);

% params
[a1, A1, C1] = get_params(y_obs, t_obs, sigma, 1);
[a2, A2, C2] = get_params(y_obs, t_obs, sigma, 2);
[a3, A3, C3] = get_params(y_obs, t_obs, sigma, 3);
[a4, A4, C4] = get_params(y_obs, t_obs, sigma, 4);
[a5, A5, C5] = get_params(y_obs, t_obs, sigma, 5);

% y_mod
y_mod1 = a1(1) + a1(2) * t_obs;
y_mod2 = a2(1) + a2(2) * t_obs + a2(3) * t_obs.^2;
y_mod3 = a3(1) + a3(2) * t_obs + a3(3) * t_obs.^2 + a3(4) * t_obs.^3;
y_mod4 = a4(1) + a4(2) * t_obs + a4(3) * t_obs.^2 + a4(4) * t_obs.^3 ...
    + a4(5) * t_obs.^4;
y_mod5 = a5(1) + a5(2) * t_obs + a5(3) * t_obs.^2 + a5(4) * t_obs.^3 ...
    + a5(5) * t_obs.^4 + a5(6) * t_obs.^5;

% dy
dy1 = ((y_obs - y_mod1) ./ sigma).^2;
dy2 = ((y_obs - y_mod2) ./ sigma).^2;
dy3 = ((y_obs - y_mod3) ./ sigma).^2;
dy4 = ((y_obs - y_mod4) ./ sigma).^2;
dy5 = ((y_obs - y_mod5) ./ sigma).^2;

% chi-square
chi1 = sum(dy1);
chi2 = sum(dy2);
chi3 = sum(dy3);
chi4 = sum(dy4);
chi5 = sum(dy5);

% y_mod plot
figure(1)
subplot(2,1,1)
hold on
plot(t_obs,y_obs,'o')
plot(t_obs,y_mod1)
plot(t_obs,y_mod2)
plot(t_obs,y_mod3)
plot(t_obs,y_mod4)
plot(t_obs,y_mod5)
title('Modellanpassung Polynome 1ter - 5ter Ordung')
xlabel('t_{i}')
ylabel('y')
legend('Messdaten', ...
    'Modellanpassung 1. Ordnung', ...
    'Modellanpassung 2. Ordnung', ...
    'Modellanpassung 3. Ordnung', ...
    'Modellanpassung 4. Ordnung', ...
    'Modellanpassung 5. Ordnung', ...
    'location','northwest')
hold off

% chi-square plot
chi = [chi1 14.067; chi2 12.592; chi3 11.070; chi4 9.488; chi5 7.815];

chi_des = {'Chi^2_1'; 'Chi^2_2'; 'Chi^2_3'; 'Chi^2_4'; 'Chi^2_5'};

subplot(2,1,2)
bar (chi)
ylim([0 40])
legend({'Chi^2_{obs}', 'Chi^2_{0.95}'});
set(gca, 'xticklabel', chi_des)
print('B1','-dpdf')

%% B2: Singulärwertzerlegung (Singular Value Decomposition) SVD
A = [1 -2  3;
    2 -1  4;
    -1 -4  1];

y  = [4 3 7]';

% Matrix regulär oder singulär
if (size(A,1) == size(A,2))
    disp('Exakt bestimmt da n == m')
    if (det(A) ~= 0)
        disp('Matrix A ist regulär da det(A) ~= 0')
    else
        disp('Matrix A ist singulär da det(A) == 0')
    end
    
    % Singulärwertzerlegung
    [U,S,V] = svd(A);
    
    % Singulärwertmatrix
    lambdaInv = inv(S);
    lambdaInv(3,3) = 0;
    
    % Berechnung der Parameter herkoemmlich vs SVD
    
    a2 = V * lambdaInv * U' * y;
    if (det(A) == 0)
        disp('Berechnung mittels inverser Matrix nicht möglich!')
        y1 = y;
    else
        a1 = A \ y;
        disp('Berechnung der Parameter mittels inverser:')
        disp(a1)
        y1 = A * a1;
    end
    
    disp('Berechnung der Parameter mittels SVD:')
    disp(a2)
    y2 = A * a2;
    error = y2 - y1;
    
    chi = sum(error.^2);
    disp('Abweichung des y-Vektors der Angabe und des berechneten y-Vektors durch Verwendung der berechneten Parameter: ')
    disp(chi)
    
    
    % Kovarianzmatrix
    lambda2 = S.^(-2);
    S2      = zeros(3);
    S2(1,1) = lambda2(1,1);
    S2(2,2) = lambda2(2,2);
    S2(3,3) = 0;
    
    Ca = V * S2 * V';
    
    disp('Kovarianzmatrix: ')
    disp(Ca)
    
elseif (size(A,1) > size(A,2))
    disp('Zeilen groesser Spalten')
    disp('Rang der Matrix: ')
    disp(rank(A))
    
elseif (size(A,1) < size(A,2))
    disp('Spalten groesser Zeilen')
    disp('Rang der Matrix: ')
    disp(rank(A))
    
end

%% B3: Parameterverteilung fuer Geradenbeispiel 1 aus Uebung E03
%Designmatrix asu Datenpunkten von B1 UE03
y_obs = [2.2; 2.6; 3.6; 3.7; 6.1; 7.6; 8.3; 8.8; 9.4];
x     = [1:9]';
sigma = 0.5;
A     = [ones(9,1),x] ./sigma;
C     = inv(A' *A);

% Parameter
sigma_d  = sqrt(C(1,1));
sigma_k  = sqrt(C(2,2));
detA     = det(A'*A);
m        = 2;
f_factor = sqrt(detA/((2*pi)^m));

delta_d = linspace(-3*sigma_d, 3*sigma_d,50);
delta_k = linspace(-3*sigma_k, 3*sigma_k,50);

% f_delta_a
f_delta_a = zeros(2,length(delta_k));
for i = 1:length(delta_d)
    for j = 1:length(delta_k)
        delta_a = [delta_d(i); delta_k(j)];
        f_delta_a(i,j) = f_factor * exp(-0.5 * delta_a' * (A' * A) * delta_a);
    end
end

% Konturlinien
con = [2.3 4.61 6.17 9.21];
cont = sqrt(detA/((2*pi)^2)).*exp(-0.5 .*con);
conh = [cont(4) cont(3) cont(2) cont(1)];

figure(2)
subplot(2,1,1)
contour(f_delta_a, conh)
title('Konturplot')
colorbar('Ticks',conh ,'TickLabels',{'68.3%','90.0%','95.4%','99.0%'})
xlabel('\Deltad')
ylabel('\Deltak')
grid on

subplot(2,1,2)
mesh(f_delta_a)
title('Parameterverteilung - 3D Plot')
xlabel('\Deltad')
ylabel('\Deltak')
zlabel('f(\Deltad,\Deltak)')
print('B3','-dpdf')