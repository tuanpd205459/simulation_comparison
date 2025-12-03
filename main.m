clc, clear, close all;

%% 1. Cài đặt tham số
addpath("D:\tuan\analysis\analysis_fringe");
%% --- 1. CẤU HÌNH VÀ KHỞI TẠO ---
num_trials = 2; % Số lần chạy
pos_coeff = cell(1, num_trials);
val_coeff = cell(1, num_trials);
noise_level = 0.3;
RMSE = zeros(num_trials, 5);
for i = 1:num_trials

    [groundtruth, pos_coeff{1,i}, val_coeff{1,i}] = create_random_zernike_surface();
    [hologram, wrapped_phase, carrier] = create_hologram(groundtruth, noise_level);
    % figure;
    % surf(groundtruth, "EdgeColor","none");
    % title("groundtruth");
    % figure;
    % imshow(hologram,[]);
    % title("anh hologram");

    estimate_phase = create_estimate_phase(hologram, groundtruth, wrapped_phase, carrier);


    % 3.3 Chạy unwrapping
    % --- Unwrap phase bằng các thuật toán ---
    phi_tie_dct      = Unwrap_TIE_DCT_Iter(wrapped_phase);     % TIE với DCT
    phi_quality      = unwrap_quality(wrapped_phase);          % Quality-guided
    phi_wls          = phase_unwrap_2dweight(wrapped_phase);      % 2D Weighted LS
    phi_goldstein    = unwrap_goldstein(wrapped_phase);     % goldstein branch-cut
    phi_proposed     = estimate_phase;                                     % Proposed / Hybrid

    value_x = size(wrapped_phase,2);
    value_y = size(wrapped_phase,1);
    phi_ls = unwrap_phase_LS(wrapped_phase, value_x, value_y);
  figure;
    surf(phi_ls, "EdgeColor","none");
    title("phi ls");
    % phi_wls = phi_ls;

    % --- Crop tất cả về cùng kích thước nhỏ nhất ---
    [groundtruth, wrapped_phase, phi_goldstein, phi_tie_dct, phi_quality, phi_wls, phi_proposed] = ...
        crop_multiple_to_smallest(groundtruth, wrapped_phase, phi_goldstein, phi_tie_dct, phi_quality, phi_wls, phi_proposed);

    phi_tie_dct = phi_tie_dct -min(phi_tie_dct(:));
    phi_quality = phi_quality -min(phi_quality(:));
    phi_wls = phi_wls -min(phi_wls(:));
    phi_goldstein = phi_goldstein -min(phi_goldstein(:));
    phi_proposed = phi_proposed -min(phi_proposed(:));
    groundtruth = groundtruth -min(groundtruth(:));

    figure;
    surf(phi_proposed, "EdgeColor","none");
    title("phi proposed");

    figure;
    surf(phi_goldstein, "EdgeColor","none");
    title("phi goldstein");

    % 3.4 Tính sai số
    RMSE(i,1) = compute_rmse(phi_goldstein, groundtruth);
    RMSE(i,2) = compute_rmse(phi_tie_dct, groundtruth);
    RMSE(i,3) = compute_rmse(phi_quality,   groundtruth);
    RMSE(i,4) = compute_rmse(phi_wls,  groundtruth);
    RMSE(i,5) = compute_rmse(phi_proposed,   groundtruth);

end
% close all;

%% 2. Thiết lập thông số vẽ
algo_names = {'Goldstein', 'TIE DCT', 'Quality', 'WLS', 'Proposed'};

% Ký hiệu marker cho Line plot để dễ phân biệt (đen trắng vẫn nhìn được)
markers = {'o', '*', 'd', '^', 's'}; 
line_styles = {'-', '-', '-', '-', '-'}; % Các đường khác nét đứt, Proposed nét liền
colors = lines(5); % Lấy bảng màu mặc định

%% 3. Vẽ biểu đồ
figure('Name', 'So sanh RMSE', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% --- BIỂU ĐỒ 1: LINE PLOT (Bên trái) ---
hold on; grid on; box on;

for i = 1:5
    if i == 5
        % Cấu hình riêng cho thuật toán Proposed (Cột 5) để nó nổi bật nhất
        plot(1:num_trials, RMSE(:, i), ...
            'LineStyle', line_styles{i}, ...
            'Marker', markers{i}, ...
            'Color', 'r', ...          % Màu đỏ
            'LineWidth', 2, ...        % Nét đậm hơn
            'MarkerFaceColor', 'b', ...
            'MarkerSize', 6);
    else
        % Các thuật toán khác
        plot(1:num_trials, RMSE(:, i), ...
            'LineStyle', line_styles{i}, ...
            'Marker', markers{i}, ...
            'Color', colors(i,:), ...
            'LineWidth', 1.5, ...
            'MarkerSize', 6);
    end
end

title('Comparison of RMSE across 15 Surfaces');
xlabel('Surface Index');
ylabel('RMSE (rad)');
xlim([1 num_trials]);
xticks(1:num_trials); % Hiển thị đủ số 1 đến 15
legend(algo_names, 'Location', 'best');
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman'); % Font chữ báo cáo

% --- BIỂU ĐỒ 2: BOX PLOT 
figure;
% Vẽ Boxplot
boxplot(RMSE, 'Labels', algo_names, 'Symbol', 'r+'); % Symbol 'r+' là dấu cộng đỏ cho outlier

title('Statistical Distribution of RMSE');
ylabel('RMSE (rad)');
grid on;
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');

% Tinh chỉnh màu sắc cho Box plot (Tùy chọn cho đẹp)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'), colors(6-j,:), 'FaceAlpha',.5);
end

%% 4. Lưu ảnh (Tùy chọn)
% exportgraphics(gcf, 'RMSE_Comparison.png', 'Resolution', 300);




function varargout = crop_multiple_to_smallest(varargin)
% Giả định tất cả các biến là 2D ma trận
n = nargin;
sizes = cellfun(@size, varargin, 'UniformOutput', false);

% Tìm kích thước nhỏ nhất theo từng chiều
min_rows = min(cellfun(@(s) s(1), sizes));
min_cols = min(cellfun(@(s) s(2), sizes));

varargout = cell(1, n);
for i = 1:n
    mat = varargin{i};
    [m, n_] = size(mat);

    % Tính chỉ số cắt đều 4 phía
    row_start = floor((m - min_rows)/2) + 1;
    col_start = floor((n_ - min_cols)/2) + 1;
    row_end = row_start + min_rows - 1;
    col_end = col_start + min_cols - 1;

    varargout{i} = mat(row_start:row_end, col_start:col_end);
end
end


