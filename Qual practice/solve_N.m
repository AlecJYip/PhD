function [N_p, N_m] = solve_N(n, phi, phi_p, pos_min)

    %eqn 1: 2 pi n N^+ - (phi +/- phi')) = pi

    %eqn 2: 2 pi n N^- - (phi +/- phi')) = -pi

    if(pos_min == 1)

        phase = phi + phi_p;

    else
        phase = phi - phi_p;

    end

    N_p = (pi+phase)/(2*pi*n);

    N_p_high = ceil(N_p);

    N_p_low = floor(N_p);

    if( abs(N_p_high-N_p) < abs(N_p-N_p_low))

        N_p = N_p_high;

    else
        N_p = N_p_low;

    end
  

    N_m = (-pi+phase)/(2*pi*n);

    N_m_high = ceil(N_m);

    N_m_low = floor(N_m);

    if( (N_m_high-N_m) < (N_m-N_m_low))

        N_m = N_m_high;

    else
        N_m = N_p_low;

    end

end