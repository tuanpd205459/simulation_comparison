function r = compute_rmse(a, b)
    r = sqrt(mean((a(:)-b(:)).^2));
end
