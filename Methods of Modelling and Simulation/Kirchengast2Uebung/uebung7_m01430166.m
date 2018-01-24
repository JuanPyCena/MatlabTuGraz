% general linear regression example - for UE Meth.Mod.Sim.
% (exercise & training example by G. Kirchengast, WEGC & IGAM/IP)
% Sonnleitner, 01430166
% 23.01.2018

clear all;
close all;
clc;

N_realizations_vec        = [10 50 300];
sigma_x                   = [0.1 0.0001];
sigma_y                   = 0.01;
sigma                     = [sigma_x(1), sigma_y];

result_X                  = [];
result_Y                  = [];
result_Xeigenvalues_check = [];
result_Aeigenvalues_check = [];
result_yestim_check       = [];
result_dyestim_check      = [];
result_Aestimate          = [];
%%
for i = 1:numel(N_realizations_vec)
    N_realizations = N_realizations_vec(i);
    xtrue = [2, 0.5, -0.25]';
    ytrue = [2.1875, 2.25, 2.1875, 2]';

    Atrue = [1, 0.5, 0.25
             1, 1  , 1
             1, 1.5, 2.25
             1, 2  , 4];

    Aeigenvalues_check = eig(Atrue' * Atrue);
    x_inv_check        = (inv(Atrue' * Atrue)*Atrue')*ytrue;
 
    X = zeros(N_realizations, 3);
    Y = zeros(N_realizations, 4);

    for n =1:N_realizations
        seed   = [n+50, n+500];
        rng(seed(1))
        dx     = sigma(1) * rand(3,1);
        rng(seed(2))
        dy     = sigma(2) * rand(4,1);
        X(n,:) = xtrue + dx;
        Y(n,:) = ytrue + dy;    
    end


    B             = (inv(X'*X)*X')*Y;
    Aestimate     = B';
    yestim_check  = Aestimate*xtrue;
    dyestim_check = yestim_check - ytrue;
    Xeigenvalues_check = eig(X'*X);
    
    result_X                  = [result_X, X];
    result_Y                  = [result_Y, Y];
    result_Xeigenvalues_check = [result_Xeigenvalues_check, Xeigenvalues_check];
    result_Aeigenvalues_check = [result_Aeigenvalues_check, Aeigenvalues_check];
    result_yestim_check       = [result_yestim_check, yestim_check];
    result_yestim_check       = [result_dyestim_check, dyestim_check];
    result_Aestimate          = [result_Aestimate, Aestimate];
    
    disp(['Rank(Aestimate): ', num2str(rank(Aestimate))])
    disp(['Rank(X): ', num2str(rank(X))])
end