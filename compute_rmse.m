function r = compute_rmse(a, b)
    diff = a - b;

    % Bỏ các giá trị NaN
    diff = diff(~isnan(diff));

    % Loại bỏ offset để so sánh công bằng
    diff = diff - mean(diff);

    % Tính RMSE
    r = sqrt(mean(diff.^2));
end
