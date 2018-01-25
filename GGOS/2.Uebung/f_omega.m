function omega = f_omega(initial, r_sun, r_moon,c_20, c_21, c_22, s_21, s_22, h, coefficient_T_g, coefficient_T_r, coefficient_F, GM_sun, GM_moon,  k_re, k_im, A, B, C, tr, timespan)

    dist_sun   = sqrt(r_sun(1,:).^2 + r_sun(2,:).^2 + r_sun(3,:).^2);
    dist_moon  = sqrt(r_moon(1,:).^2 + r_moon(2,:).^2 + r_moon(3,:).^2);

    M_sun(1,:) = (((C-B) .* r_sun(2,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(2,:) = (((A-C) .* r_sun(1,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(3,:) = (((B-A) .* r_sun(1,:) .* r_sun(2,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);

    M_moon(1,:) = (((C-B) .* r_moon(2,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(2,:) = (((A-C) .* r_moon(1,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(3,:) = (((B-A) .* r_moon(1,:) .* r_moon(2,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);

    M = M_sun + M_moon;

    T_g_c1 = coefficient_T_g .* [((1/sqrt(3)) .* c_20 - c_22); -s_22; -c_21];
    T_g_c2 = coefficient_T_g .* [-s_22; ((1/sqrt(3)) .* c_20 + c_22); -s_21];
    T_g_c3 = coefficient_T_g .* [-c_21; -s_21; -(2/sqrt(3)) .* c_20];

    T_g_c1(1,:) = T_g_c1(1,:) + tr/3;
    T_g_c2(2,:) = T_g_c2(2,:) + tr/3;
    T_g_c3(3,:) = T_g_c3(3,:) + tr/3;

    T_g_c1(:,end+1) = T_g_c1(:,end);
    T_g_c2(:,end+1) = T_g_c2(:,end);
    T_g_c3(:,end+1) = T_g_c3(:,end);

    h(:,end+1) = h(:,end);

    for time = 1:timespan+1
        T_g{time}  = [T_g_c1(:, time),T_g_c2(:, time), T_g_c3(:, time)];
    end

    for time = 1:timespan
        dT_g{time} = T_g{time + 1} - T_g{time};
    end

    for time = 1:(length(h)-1)
        dh{time} = h(:,time + 1) - h(:,time);
    end


    %% Calculation
    dt = 1;
    counter = 1;
    w = initial;
    for t = 1:dt:24*365

        w_x = w(1);
        w_y = w(2);
        w_z = w(3);

        T_r = coefficient_T_r .* [0,     0,      k_re * w_x + k_im * w_y;
                                  0,     0,      k_re * w_y - k_im * w_x;
                k_re * w_x + k_im * w_y,     k_re * w_y - k_im * w_x,         0];

        F = T_g{t} + T_r + coefficient_F .* [k_re, k_im, 0;
                                            -k_im, k_re, 0;
                                             0   ,  0  , 0];

        w = w + dt .* (F\(M(:,t) - (dT_g{t} * w) - cross(w, (T_g{t} + T_r) * w) - cross(w,h(:,t)) - dh{t})); 

        result{t} = w;
    end

    omega = cell2mat(result);

end

