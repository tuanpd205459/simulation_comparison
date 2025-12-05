clc, clear, close all;

%% 1. Cài đặt tham số
addpath("D:\tuan\analysis\analysis_fringe");
%% --- 1. CẤU HÌNH VÀ KHỞI TẠO ---
%% 1. Thiết lập chung và chạy vòng lặp
num_trials = 10; % Số lần chạy
RMSE = zeros(num_trials, 4); % 4 thuật toán: Goldstein, Quality, WLS, Proposed
noise_level = 0.08;
% Cấu hình tên và màu sắc chuẩn bị cho vẽ
algo_names = {'Goldstein', 'Quality', 'WLS', 'Proposed'};

% --- Cấu hình Style cho Publication ---
% Màu sắc: [R G B]
c_gold   = [0 0 1];       % Blue
c_qual   = [0 0.5 0];     % Dark Green
c_wls    = [0.2 0.2 0.2]; % Dark Gray (Black-ish)
c_prop   = [0 1 1];       % CYAN (Theo yêu cầu)

colors = [c_gold; c_qual; c_wls; c_prop];

% Markers và Line Styles
markers     = {'^', 's', 'd', 'o'}; 
line_styles = {'--', '--', '--', '-'}; % Proposed nét liền, còn lại nét đứt
line_widths = [1.2, 1.2, 1.2, 2.0];    % Proposed dày hơn

for i = 1:num_trials
    % --- 1. Tạo dữ liệu mô phỏng ---
    [groundtruth, ~, ~] = create_random_zernike_surface();
    [hologram, wrapped_phase, carrier, object_with_noise] = create_hologram(groundtruth, noise_level);
    estimate_phase = create_estimate_phase(hologram, groundtruth, wrapped_phase, carrier);

    % --- 2. Chạy thuật toán Unwrapping ---
    phi_goldstein = unwrap_goldstein(wrapped_phase);
    phi_quality   = unwrap_quality(wrapped_phase);
    phi_wls       = phase_unwrap_2dweight(wrapped_phase);
    phi_proposed  = estimate_phase; 

    % (Bỏ qua phi_ls vì không thấy trong danh sách algo_names)

    % --- 3. Crop về cùng kích thước ---
    [groundtruth, object_with_noise, wrapped_phase, phi_goldstein, phi_quality, phi_wls, phi_proposed] = ...
        crop_multiple_to_smallest(groundtruth, object_with_noise, wrapped_phase, phi_goldstein, phi_quality, phi_wls, phi_proposed);

    % --- 4. Xử lý Offset (Quan trọng: Trừ Mean để tính RMSE chính xác) ---
    % Lưu ý: Khi tính RMSE cho bài báo, ta thường so sánh dao động quanh giá trị trung bình
    % hoặc phải shift sao cho phase khớp nhau nhất. Ở đây dùng trừ mean.
    gt_zero   = object_with_noise - mean(object_with_noise(:));
    p_gold_z  = phi_goldstein - mean(phi_goldstein(:));
    p_qual_z  = phi_quality   - mean(phi_quality(:));
    p_wls_z   = phi_wls       - mean(phi_wls(:));
    p_prop_z  = phi_proposed  - mean(phi_proposed(:));

    % --- 5. Tính RMSE và lưu vào ma trận (Cột 1->4) ---
    RMSE(i, 1) = compute_rmse(p_gold_z, gt_zero); % Goldstein
    RMSE(i, 2) = compute_rmse(p_qual_z, gt_zero); % Quality
    RMSE(i, 3) = compute_rmse(p_wls_z,  gt_zero); % WLS
    RMSE(i, 4) = compute_rmse(p_prop_z, gt_zero); % Proposed
end

%% 2. Vẽ hình (Visualization)
% Tạo một Figure lớn chứa cả 2 biểu đồ (hoặc tách ra nếu muốn)
figure('Name', 'Scientific RMSE Comparison', 'Color', 'w', 'Position', [50, 50, 1000, 500]);

% --- BIỂU ĐỒ 1: LINE PLOT (Theo dõi qua các lần chạy) ---
subplot(1, 2, 1);
hold on; grid on; box on;

for k = 1:4
    if k == 4 % Cấu hình riêng cho Proposed (Cyan)
        plot(1:num_trials, RMSE(:, k), ...
            'LineStyle', line_styles{k}, ...
            'Marker', markers{k}, ...
            'Color', colors(k, :), ...          % Màu Cyan
            'LineWidth', line_widths(k), ...
            'MarkerFaceColor', 'w', ...         % Mặt trong marker màu trắng cho nổi
            'MarkerEdgeColor', [0 0.5 0.5], ... % Viền marker tối hơn xíu để rõ trên nền Cyan
            'MarkerSize', 7);
    else
        plot(1:num_trials, RMSE(:, k), ...
            'LineStyle', line_styles{k}, ...
            'Marker', markers{k}, ...
            'Color', colors(k, :), ...
            'LineWidth', line_widths(k), ...
            'MarkerSize', 6);
    end
end

% Trang trí biểu đồ 1
title('RMSE Performance per Trial');
xlabel('Surface Index (Trial)');
ylabel('RMSE (rad)');
xlim([0.5, num_trials + 0.5]);
xticks(1:max(1, round(num_trials/5)):num_trials); % Hiển thị tick trục X thông minh
legend(algo_names, 'Location', 'best', 'FontSize', 10);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);


% --- BIỂU ĐỒ 2: BOX PLOT (Phân bố thống kê) ---
subplot(1, 2, 2);
box_handle = boxplot(RMSE, 'Labels', algo_names, 'Symbol', 'rx'); % 'rx' là outlier dấu x đỏ

% Trang trí Boxplot
title('Statistical Distribution of RMSE');
ylabel('RMSE (rad)');
grid on; box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);

% --- Tô màu Boxplot cho đồng bộ với Lineplot ---
h = findobj(gca, 'Tag', 'Box');
% Boxplot vẽ ngược từ dưới lên, nên ta duyệt ngược
for j = 1:length(h)
    algo_idx = length(h) - j + 1; % Map ngược lại index thuật toán
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(algo_idx, :), ...
          'FaceAlpha', 0.3, 'EdgeColor', colors(algo_idx, :), 'LineWidth', 1.2);
end

% Thêm tiêu đề chung (nếu dùng MATLAB bản mới)
% sgtitle('Performance Comparison of Phase Unwrapping Algorithms', 'FontName', 'Times New Roman', 'FontWeight', 'bold');

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


