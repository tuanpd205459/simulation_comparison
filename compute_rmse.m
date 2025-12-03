function r = compute_rmse(a, b)
    diff = a - b;
    r = sqrt(mean(diff(~isnan(diff)).^2));
end
