function [T_g, w] = calculate_w(x, dt, w_in, T_g0, h, r_sun, r_moon, dh, coefficient_F,coefficient_T_g, coefficient_T_r, GM_sun, GM_moon, delta_x)
    c_20 = x(1);
    c_21 = x(2);
    c_22 = x(3);
    s_21 = x(4);
    s_22 = x(5);
    k_re = x(6);
    k_im = x(7);
    A    = x(8);
    B    = x(9);
    C    = x(10);

    delta_c20 = delta_x(1);
    delta_c21 = delta_x(2);
    delta_c22 = delta_x(3);
    delta_s21 = delta_x(4);
    delta_s22 = delta_x(5);
    delta_kre = delta_x(6);
    delta_kim = delta_x(7);
    delta_A   = delta_x(8);
    delta_B   = delta_x(9);
    delta_C   = delta_x(10);
    
    c_20 = c_20 + delta_c20;
    c_21 = c_21 + delta_c21;
    c_22 = c_22 + delta_c22;
    s_21 = s_21 + delta_s21;
    s_22 = s_22 + delta_s22;
    k_re = k_re + delta_kre;
    k_im = k_im + delta_kim;
    A    = A    + delta_A;
    B    = B    + delta_B;
    C    = C    + delta_C;
    
    w_x = w_in(1);
    w_y = w_in(2);
    w_z = w_in(3);
    
    tr = A + B + C;
    
    dist_sun   = sqrt(r_sun(1).^2 + r_sun(2).^2 + r_sun(3).^2);
    dist_moon  = sqrt(r_moon(1).^2 + r_moon(2).^2 + r_moon(3).^2);

    M_sun(1,:) = (((C-B) .* r_sun(2) .* r_sun(3)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(2,:) = (((A-C) .* r_sun(1) .* r_sun(3)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(3,:) = (((B-A) .* r_sun(1) .* r_sun(2)) ./ (dist_sun.^5)) .* (3 * GM_sun);

    M_moon(1,:) = (((C-B) .* r_moon(2) .* r_moon(3)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(2,:) = (((A-C) .* r_moon(1) .* r_moon(3)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(3,:) = (((B-A) .* r_moon(1) .* r_moon(2)) ./ (dist_moon.^5)) .* (3 * GM_moon);

    M = M_sun + M_moon;
    
    T_g_c1 = coefficient_T_g .* [((1/sqrt(3)) .* c_20 - c_22); -s_22; -c_21];
    T_g_c2 = coefficient_T_g .* [-s_22; ((1/sqrt(3)) .* c_20 + c_22); -s_21];
    T_g_c3 = coefficient_T_g .* [-c_21; -s_21; -(2/sqrt(3)) .* c_20];
    
    T_g_c1 = T_g_c1(1) + tr/3;
    T_g_c2 = T_g_c2(2) + tr/3;
    T_g_c3 = T_g_c3(3) + tr/3;
    
    T_g = [T_g_c1, T_g_c2, T_g_c3];
    
    dT_g = T_g - T_g0;
    
    T_r = coefficient_T_r .* [0,     0,      k_re * w_x + k_im * w_y;
                              0,     0,      k_re * w_y - k_im * w_x;
            k_re * w_x + k_im * w_y,     k_re * w_y - k_im * w_x,         0];
        
    F = T_g + T_r + coefficient_F .* [k_re, k_im, 0;
                                     -k_im, k_re, 0;
                                     0   ,  0  , 0];
    
    w = w_in + dt .* (F\(M - (dT_g * w_in) - cross(w_in, (T_g + T_r) * w_in) - cross(w_in,h) - dh)); 
end

