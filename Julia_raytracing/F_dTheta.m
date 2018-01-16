function dTheta= F_dTheta(phi0,r,step)
%Parameter
N0 = 315;
Hn = 7;
r0 = 6371;



Nr = N0 * exp(-((r+step) - r0) / Hn);         %Refraktivität
nr = 1 + 10e-6 * Nr;                    %Brechungsindexfeld d. Atmosphäre
C = 1 * r0 * sind(phi0);             %Konstante
phir = asind(C ./ (nr .* (r+step)));

dTheta = C ./ (nr .* (r+step).^2 .* cosd(phir));
end

