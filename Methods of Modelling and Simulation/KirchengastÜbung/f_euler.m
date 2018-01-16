function [ret_theta, ret_r] = f_euler(phi0, Hn, N0, theta, begin, step, stop)
    r0          = begin;
    theta_array = [theta];
    r_array     = [r0];
    
    for r = begin:step:stop         
        N_r   = N0 * exp(-(r-r0)/Hn);
        n_r   = 1 + 1e-6 * N_r;
        C     = (1+1e-6*N0) * r0 * sind(phi0);
        phi_r = asind(C/(n_r * r));

        dtheta = C/(n_r * r * r * cos(phi_r));

        k = step *dtheta;
        theta = theta + k;
        
        theta_array = [theta_array, theta];
        r_array     = [r_array, r];
    end
    
    ret_theta = theta_array;
    ret_r     = r_array;
end

