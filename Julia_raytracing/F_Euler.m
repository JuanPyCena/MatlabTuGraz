function [theta,r] = F_Euler(phi0)
%Euler Berechnet den Strahlengang in der Atmosph�re nach der Eulerschen
%Methode, abh�ngig von Phi

%Parameter
r0 = 6371;
N0 = 315;
Hn = 7;
r = r0 : 0.1 : r0+650;
theta(1) = 0;

Nr = N0 * exp(-(r - r0) / Hn);         %Refraktivit�t
nr = 1 + 10e-6 * Nr;                    %Brechungsindexfeld d. Atmosph�re
C = 1 * r0 * sind(phi0);             %Konstante
phir = asind(C ./ (nr .* r));

%Euler
dTheta = C ./ (nr .* r.^2 .* cosd(phir));
for i = 1:length(dTheta)
   theta(i+1) = theta(i) + 0.1 * dTheta(i); 
end


end

