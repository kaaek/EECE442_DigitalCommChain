function E = errorEnergy (xt, xhat, t)
    integrand = (xt - xhat).^2;
    E = integrate(t, integrand);
end