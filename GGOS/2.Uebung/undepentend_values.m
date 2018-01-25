function [initial,G, GM_sun, GM_moon, omega_N, Mass, R, A, B, C, coefficient_F, coefficient_T_g, coefficient_T_r, M] = undepentend_values(r_moon, r_sun, reference)
    
    initial = reference(:,1) .*3600;              % [rad / h]
    G       = (6.674e-11) * (3600 * 3600);        % [m^3/(kg* h^2)]
    GM_sun  = (1.32712442076e20) * (3600 * 3600); % [m^3/ h^2]
    GM_moon = (4.9027779e12) * (3600 * 3600);     % [m^3/ h^2]
    omega_N = (7.2921151467064e-5) * 3600;        % [rad/h]
    Mass    = 5.9737e24;                          % [kg]
    R       = 6378136.6;                          % [m]
    A       = 0.3296108 * Mass * R * R;           % [kg * m^2]
    B       = 0.3296108 * Mass * R * R;           % [kg * m^2]
    C       = 0.3307007 * Mass * R * R;           % [kg * m^2]

    coefficient_T_g = sqrt(5/3) * Mass * R * R;
    coefficient_T_r = (omega_N * R^5) / (3 * G);
    coefficient_F   = (omega_N * omega_N * R^5) / (3 * G);

    dist_sun   = sqrt(r_sun(1,:).^2 + r_sun(2,:).^2 + r_sun(3,:).^2);
    dist_moon   = sqrt(r_moon(1,:).^2 + r_moon(2,:).^2 + r_moon(3,:).^2);

    M_sun(1,:) = (((C-B) .* r_sun(2,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(2,:) = (((A-C) .* r_sun(1,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(3,:) = (((B-A) .* r_sun(1,:) .* r_sun(2,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);

    M_moon(1,:) = (((C-B) .* r_moon(2,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(2,:) = (((A-C) .* r_moon(1,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(3,:) = (((B-A) .* r_moon(1,:) .* r_moon(2,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    
    M = M_sun + M_moon;
end

