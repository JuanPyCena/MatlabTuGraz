function [ret_theta, ret_r] = f_euler(phi0, Hn, N0, theta, begin, step, stop)
    r0          = begin;
    theta_array = theta;
    r_array     = r0;
    
    for r = begin:step:stop 
        
        k = step * f_dtheta(r, r0, phi0, Hn, N0);
        theta = theta + k;
        
        theta_array = [theta_array, theta];
        r_array     = [r_array, r];
    end
    
    ret_theta = theta_array;
    ret_r     = r_array;
end

