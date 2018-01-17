function dtheta = f_dtheta(r, r0, phi0, Hn, N0)
    N_r   = N0 * exp(-(r-r0)/Hn);
    n_r   = 1 + 1e-6 * N_r;
    C     = (1+1e-6*N0) * r0 * sin(phi0);
    phi_r = asin(C/(n_r * r));

    dtheta = C/(n_r * r * r * cos(phi_r));
end

