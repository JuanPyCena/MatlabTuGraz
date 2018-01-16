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
a = 1
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

N_ex = explicit(Ni0,1);
figure(1)
mesh(1:100:100*length(N_ex),100:600,N_ex ./N00)
title('Explizite Lösungsform, \Deltat_{krit}')
xlabel('Zeit [s]')
ylabel('Höhe [km]')
zlabel('n_e / N_{00}')

N_ex02 = explicit(Ni0,0.2);
figure(2)
mesh(1:100:100*length(N_ex02),100:600,N_ex02 ./N00)
title('Explizite Lösungsform, 0.2\Deltat_{krit}')
xlabel('Zeit [s]')
ylabel('Höhe [km]')
zlabel('n_e / N_{00}')

N_ex2 = explicit(Ni0,2);
figure(3)
plot(100:600,N_ex2(:,1) ./N00)
title('Explizite Lösungsform, 2\Deltat_{krit}')
xlabel('Höhe [m]')
ylabel('n_e / N_{00}')
hold on
plot(100:600,N_ex2(:,5) ./N00,'c')
plot(100:600,N_ex2(:,7) ./N00,'m')
hold off
legend('t = 1','t = 5','t = 7')

%% implizit
n = zeros(size(Ni0,1),tmax/dt);
for j = 1:length(Ni0)
    n(j,1) = Ni0(j);
end
alpha = (A0 * a * dt)/((db*1000)^2).*ones(size(n,1),1)';
alpha1 = (A0 * a * dt)/((db*1000)^2).*ones(size(n,1)-1,1)';

B = diag(-alpha1,1);
C = diag(-alpha1,-1);
A = diag(1 + 2*alpha)+ B + C;
A(1,1) = 1 - delta0;
A(1,2) = 0;
A(end,end-1) = 0;
A(end,end) = 1-delta1;


for t = 1:tmax/dt
    n(:,t + 1) = A\n(:,t);
    i_max = 251;
    index_a = find(n(1:i_max,t+1) < N00);
    n(index_a,t+1) = N00;
    index_e = find(n(i_max:end,t+1) < N10);
    n(index_e + i_max-1,t) = N10;
end

mesh(1:100:length(n)*100,100:600,n ./N00)
title('Implizite Lösungsform, 2 \Deltat_{krit}')
xlabel('Zeit [s]')
ylabel('Höhe [km]')
zlabel('n_e / N_{00}')


diff = n - N_ex02;

mesh(1:100:length(n)*100,100:600,diff./N00)
title('Differenz Implizit explizit, \Deltat_{krit}')
xlabel('Zeit [s]')
ylabel('Höhe [km]')
zlabel('n_e / N_{00}')