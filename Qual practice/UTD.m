function [D_par, D_per] = UTD(L, phi, phi_p, n, beta,  gamma_p)

    term_1 = -exp(-1i*pi/4)...
        /(2*n*sqrt(2*pi*beta)*sin(gamma_p));

    pos_min = -1;

    [N_p_m_1, N_m_m_1] = solve_N(n, phi, phi_p, pos_min);

    [a_p_m_1, a_m_m_1] = solve_a(n, phi, phi_p, pos_min, N_p_m_1, N_m_m_1);


    pos_min = 1;

    [N_p_1, N_m_1] = solve_N(n, phi, phi_p, pos_min);

    [a_p_1, a_m_1] = solve_a(n, phi, phi_p, pos_min, N_p_1, N_m_1);



    term_2 = cot((pi+(phi-phi_p))/(2*n))*fresnel(beta, L, a_p_m_1);

    term_3 = cot((pi-(phi-phi_p))/(2*n))*fresnel(beta, L, a_m_1);



    term_4 = cot((pi+(phi+phi_p))/(2*n))*fresnel(beta, L, a_p_1);

    term_5 = cot((pi-(phi+phi_p))/(2*n))*fresnel(beta, L, a_m_m_1);

    D_par = term_1*(term_2+term_3-term_4-term_5);

    D_per = term_1*(term_2+term_3+term_4+term_5);

end