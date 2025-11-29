function foo = MSE (xt, xhat, t)
    foo = errorEnergy(xt, xhat, t)/length(t);
end