function [n] = indexOfRefrection(r)
N0 = 315;
r0 = 6371; %[km]
Hn = 7; %Skalenhohe [km]
N = N0.*exp(-(r-r0)./Hn);
n = 1 + (1E-6)*N;