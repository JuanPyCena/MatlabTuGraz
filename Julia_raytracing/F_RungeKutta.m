function [theta,r] = F_RungeKutta(phi0)
%RungeKutta Berechnet den Strahlengang in der Atmosphäre nach der
%Runge-Kutta Methode, abhängig von Phi

dh = 0.1;
theta(1) = 0;
r = 6371 : 0.1 : 6371+650;
k1 = dh * F_dTheta(phi0, r,0);
k2 = dh * F_dTheta(phi0, r, dh/2);
k3 = dh * F_dTheta(phi0, r, dh/2);
k4 = dh * F_dTheta(phi0, r, dh);

for i = 1:length(k1)
   theta(i+1) = theta(i) + k1(i)/6 + k2(i)/3 + k3(i)/3 + k4(i)/6;
end

end

