function [ret_theta, ret_r] = f_rungeKutta(phi0, Hn, N0, theta, begin, step, stop)
    r0          = begin;
    theta_array = theta;
    r_array     = r0;
    
    for r = begin:step:stop          
        k1 = step * f_dtheta(r         , r0, phi0, Hn, N0);
        k2 = step * f_dtheta(r + step/2, r0, phi0, Hn, N0);
        k3 = step * f_dtheta(r + step/2, r0, phi0, Hn, N0);
        k4 = step * f_dtheta(r + step  , r0, phi0, Hn, N0);
        
        theta = theta + k1/6 + k2/3 + k3/3 + k4/6;
        
        theta_array = [theta_array, theta];
        r_array     = [r_array, r];
    end
    
    ret_theta = theta_array;
    ret_r     = r_array;
end

