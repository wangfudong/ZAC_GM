function sLen = AdaptiveStep(x, x_est, sLen, theta, rho)
    N = sqrt(size(x, 1));
    if trace((x_est - x)' * (x_est - x)) < theta * N
        sLen = rho * sLen;
    else
        sLen = max(1, sLen / rho);
    end
end