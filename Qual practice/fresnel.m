function F = fresnel(beta, L, a_pm)

    X = beta*L*a_pm;

    int_start = abs(sqrt(X));

    int_end = 5000;

    d_X = 0.0001;

    tau = int_start:d_X:int_end;

    F = 2*1i*abs(sqrt(X))*exp(1i*X)*sum(exp(-1i *tau.^2)*d_X); % Vectorized summation

end