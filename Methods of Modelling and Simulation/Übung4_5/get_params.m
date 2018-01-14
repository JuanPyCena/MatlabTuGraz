function [a, A, C] = get_params(y_obs, x_obs, sigma, order)
% Uebung E04/E05 Method of Modelling and Simulation
% Sonnleitner Felix, 01430166
% 09.12.2017

% Initalizing
a = zeros(order,1);
A = [ones(length(y_obs),1) ./sigma];

counter = 1;

% Creating Design matrix
while (counter <= order)  
    mod = (x_obs.^counter)./sigma;
    A   = [A, mod];
    counter = counter + 1;
end

C = inv(A' * A);

a = inv(A' * A) * A' * (y_obs./sigma);
end

