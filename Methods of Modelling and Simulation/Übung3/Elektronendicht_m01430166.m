clear all;
clc;
close all;

A0 = 5e3;
N00 = 2e10;
N10 = 20e10;
delta0 = 1e-5 * N00;
delta1 = 1e-5 * N10;

b0 = 100;
b1 = 600;
db = 1;
b = b0:db:b1;

t0 = 0;
tmax = 172800;
dt = (db*1000)^2/(2 * A0)
t = t0:dt:tmax;

Nmax = 10e11;
bmax = 350;
wN = 10;
Ni0 = Nmax * exp(-((b-bmax).^2)/wN^2);
Ni0 = Ni0';

for i = 1:length(Ni0)
    if (i < length(Ni0)/2)
        if (Ni0(i) < N00)
            Ni0(i) = N00;
        end
    end
    if (i >= (length(Ni0)/2))
        if (Ni0(i) < N10)
            Ni0(i) = N10;
        end
    end
end
plot(Ni0)
% 
% %% explizit
% format long
% n = zeros(length(Ni0),tmax/dt);
% for j = 1:length(Ni0)
%     n(j,1) = Ni0(j);
% end
% 
% for t = 1:30
%     t
%     for k = 2:501
%         
%         k
%         n(k,t+1) = n(k,t) + dt*A0(n(k+1,t) - 2 * n(k,t) + n(k-1,t))  
%     end
% end

%% implizit
n = zeros(size(Ni0,1),tmax/dt);
for j = 1:length(Ni0)
    n(j,1) = Ni0(j);
end
alpha = (A0 * dt)/((db*1000)^2).*ones(size(n,1)-1,1)';
B = diag(-alpha,1);
C = diag(-alpha,1);
A = (1 + 2*alpha(1)) .* (ones(size(n,1))) + B + C;
A(1,1) = 1 - delta0;
A(end,end) = 1-delta1;

for t = 1:tmax/dt
    n(:,t + 1) = A\n(:,t);
end

%% plot

mesh(n)