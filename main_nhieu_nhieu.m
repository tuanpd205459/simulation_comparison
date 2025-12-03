clc, clear, close all;
%% 1. Cài đặt đường dẫn và tham số
addpath("D:\tuan\analysis\analysis_fringe");

%% --- CẤU HÌNH ---
% 1. Tạo MỘT bề mặt gốc (Ground Truth) duy nhất bên ngoài vòng lặp
% Nếu hàm create_random_zernike_surface không có sẵn, có thể thay bằng peaks(256)
[groundtruth, pos_coeff, val_coeff] = create_random_zernike_surface();

figure;
surf(groundtruth, "EdgeColor", "none");
title("Ground Truth Gốc");
axis tight;

% 2. Thiết lập dải nhiễu muốn kiểm tra
noise_level_range = 0.08 : 0.02 : 0.15; 
num_trials = length(noise_level_range);

% Khởi tạo ma trận lưu kết quả RMSE
% Cột: 1-Goldstein, 2-TIE DCT, 3-Quality, 4-WLS, 5-Proposed
RMSE = zeros(num_trials, 5); 

fprintf('Dang chay mo phong voi 1 be mat, thay doi %d muc nhieu...\n', num_trials);

%% --- VÒNG LẶP THEO MỨC ĐỘ NHIỄU ---
for i = 1:num_trials
    fprintf("Dang chay lan: %d / %d (Noise: %.2f)\n", i, num_trials, noise_level_range(i));
    
    current_noise = noise_level_range(i);
    
    [hologram, wrapped_phase, carrier] = create_hologram(groundtruth, current_noise);

    estimate_phase = create_estimate_phase(hologram, groundtruth, wrapped_phase, carrier);
    
    % (Tùy chọn) Chỉ vẽ hình ở lần chạy đầu tiên hoặc cuối cùng để tránh quá nhiều cửa sổ
    if i == 1 
        figure; imshow(hologram, []); title(['Hologram (Noise ' num2str(current_noise) ')']);
        figure; surf(carrier, "EdgeColor","none"); title('Anh carrier');
        fft_img = fftshift(fft2(hologram));
        figure; imshow(log(abs(fft_img) + 1), []); title('Pho bien do FFT');
    end

    % --- 3.3 Chạy các thuật toán Unwrapping ---
    phi_goldstein = unwrap_goldstein(wrapped_phase);
    phi_tie_dct   = Unwrap_TIE_DCT_Iter(wrapped_phase);
    phi_quality   = unwrap_quality(wrapped_phase);
    phi_wls       = phase_unwrap_2dweight(wrapped_phase);
    phi_proposed  = estimate_phase;

    % --- Xử lý cắt ảnh (Crop) để cùng kích thước ---
    % QUAN TRỌNG: Đổi tên biến đầu ra thành *_cut để KHÔNG ghi đè lên biến groundtruth gốc
    [groundtruth_cut, wrapped_phase_cut, phi_goldstein_cut, phi_tie_cut, phi_qual_cut, phi_wls_cut, phi_proposed_cut] = ...
        crop_multiple_to_smallest(groundtruth, wrapped_phase, phi_goldstein, phi_tie_dct, phi_quality, phi_wls, phi_proposed);

    % Trừ trung bình (DC offset) cho các biến ĐÃ CẮT
    groundtruth_cut   = groundtruth_cut - mean(groundtruth_cut(:));
    phi_goldstein_cut = phi_goldstein_cut - mean(phi_goldstein_cut(:));
    phi_tie_cut       = phi_tie_cut - mean(phi_tie_cut(:));
    phi_qual_cut      = phi_qual_cut - mean(phi_qual_cut(:));
    phi_wls_cut       = phi_wls_cut - mean(phi_wls_cut(:));
    phi_proposed_cut  = phi_proposed_cut - mean(phi_proposed_cut(:));

    % --- 3.4 Tính sai số RMSE trên các biến ĐÃ CẮT ---
    RMSE(i,1) = compute_rmse(phi_goldstein_cut, groundtruth_cut);
    RMSE(i,2) = compute_rmse(phi_tie_cut,       groundtruth_cut);
    RMSE(i,3) = compute_rmse(phi_qual_cut,      groundtruth_cut);
    RMSE(i,4) = compute_rmse(phi_wls_cut,       groundtruth_cut);
    RMSE(i,5) = compute_rmse(phi_proposed_cut,  groundtruth_cut);
end

%% 2. Thiết lập thông số vẽ
algo_names = {'Goldstein', 'TIE DCT', 'Quality', 'WLS', 'Proposed'};
markers = {'o', '*', 'd', '^', 's'}; 
line_styles = {'--', '--', '--', '--', '-'}; % Proposed nét liền
colors = lines(5); 

%% 3. Vẽ biểu đồ so sánh
figure('Name', 'RMSE vs Noise Level', 'Color', 'w', 'Position', [100, 100, 1000, 600]);

% --- BIỂU ĐỒ LINE PLOT ---
hold on; grid on; box on;

for k = 1:5
    if k == 5
        % Cấu hình nổi bật cho Proposed
        plot(noise_level_range, RMSE(:, k), ...
            'LineStyle', line_styles{k}, ...
            'Marker', markers{k}, ...
            'Color', 'r', 'LineWidth', 2.5, ...
            'MarkerFaceColor', 'r', 'MarkerSize', 8);
    else
        % Các thuật toán khác
        plot(noise_level_range, RMSE(:, k), ...
            'LineStyle', line_styles{k}, ...
            'Marker', markers{k}, ...
            'Color', colors(k,:), 'LineWidth', 1.5, ...
            'MarkerSize', 6);
    end
end

title('Performance under varying Noise Levels (Single Surface)');
xlabel('Noise Level (Standard Deviation)');
ylabel('RMSE (rad)');
legend(algo_names, 'Location', 'northwest');
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
xlim([min(noise_level_range), max(noise_level_range)]);


% %% 1. Cài đặt đường dẫn và tham số
% addpath("D:\tuan\analysis\analysis_fringe");
% 
% %% --- CẤU HÌNH ---
% % 1. Tạo MỘT bề mặt gốc (Ground Truth) duy nhất bên ngoài vòng lặp
% [groundtruth, pos_coeff, val_coeff] = create_random_zernike_surface();
% figure;
% surf(groundtruth,"EdgeColor","none");
% title("ground truth");
% % 2. Thiết lập dải nhiễu muốn kiểm tra (Ví dụ: từ 0.1 đến 0.5)
% noise_level_range = 0.08 : 0.01 : 0.11; 
% num_trials = length(noise_level_range);
% 
% % Khởi tạo ma trận lưu kết quả RMSE
% % Cột: 1-Goldstein, 2-TIE DCT, 3-Quality, 4-WLS, 5-Proposed
% RMSE = zeros(num_trials, 5); 
% 
% fprintf('Dang chay mo phong voi 1 be mat, thay doi %d muc nhieu...\n', num_trials);
% 
% %% --- VÒNG LẶP THEO MỨC ĐỘ NHIỄU ---
% for i = 1:num_trials
%     fprintf("Dang chay lan: %d \n", i);
%     current_noise = noise_level_range(i);
%     [hologram, wrapped_phase, carrier] = create_hologram(groundtruth, current_noise);
% 
%     estimate_phase = create_estimate_phase(hologram, groundtruth, wrapped_phase, carrier);
%     figure;
%     imshow(hologram,[]);
%     title("anh hologram");
% 
%     figure;
%     surf(carrier, "EdgeColor","none");
%     title("anh carrier");
%     fft = fftshift(fft2(hologram));
%     figure;
%     imshow(log(abs(fft) + 1), []);
%     title("Pho bien do FFT");
% 
% 
%     % --- 3.3 Chạy các thuật toán Unwrapping ---
%     phi_goldstein = unwrap_goldstein(wrapped_phase);
%     phi_tie_dct = Unwrap_TIE_DCT_Iter(wrapped_phase);
%     phi_quality = unwrap_quality(wrapped_phase);
%     phi_wls = phase_unwrap_2dweight(wrapped_phase);
%     phi_proposed = estimate_phase;
% 
%     % --- Xử lý cắt ảnh (Crop) để cùng kích thước ---
%     [groundtruth, wrapped_phase, phi_goldstein, phi_tie, phi_qual, phi_wls_c, phi_proposed] = ...
%         crop_multiple_to_smallest(groundtruth, wrapped_phase, phi_goldstein, phi_tie_dct, phi_quality, phi_wls, phi_proposed);
% 
%     groundtruth   = groundtruth - mean(groundtruth(:));
%     phi_goldstein  = phi_goldstein - mean(phi_goldstein(:));
%     phi_tie   = phi_tie - mean(phi_tie(:));
%     phi_qual  = phi_qual - mean(phi_qual(:));
%     phi_wls_c = phi_wls_c - mean(phi_wls_c(:));
%     phi_proposed  = phi_proposed - mean(phi_proposed(:));
% 
%     % --- 3.4 Tính sai số RMSE ---
%     RMSE(i,1) = compute_rmse(phi_goldstein,  groundtruth);
%     RMSE(i,2) = compute_rmse(phi_tie,   groundtruth);
%     RMSE(i,3) = compute_rmse(phi_qual,  groundtruth);
%     RMSE(i,4) = compute_rmse(phi_wls_c, groundtruth);
%     RMSE(i,5) = compute_rmse(phi_proposed,  groundtruth);
% end
% 
% %% 2. Thiết lập thông số vẽ
% algo_names = {'Goldstein', 'TIE DCT', 'Quality', 'WLS', 'Proposed'};
% markers = {'o', '*', 'd', '^', 's'}; 
% line_styles = {'--', '--', '--', '--', '-'}; % Proposed nét liền, còn lại nét đứt
% colors = lines(5); 
% 
% %% 3. Vẽ biểu đồ so sánh
% figure('Name', 'RMSE vs Noise Level', 'Color', 'w', 'Position', [100, 100, 1000, 600]);
% 
% % --- BIỂU ĐỒ 1: LINE PLOT (RMSE theo mức nhiễu) ---
% hold on; grid on; box on;
% 
% for k = 1:5
%     if k == 5
%         % Cấu hình nổi bật cho Proposed
%         plot(noise_level_range, RMSE(:, k), ... % Trục X là noise_level_range
%             'LineStyle', line_styles{k}, ...
%             'Marker', markers{k}, ...
%             'Color', 'r', 'LineWidth', 2.5, ...
%             'MarkerFaceColor', 'r', 'MarkerSize', 8);
%     else
%         % Các thuật toán khác
%         plot(noise_level_range, RMSE(:, k), ... % Trục X là noise_level_range
%             'LineStyle', line_styles{k}, ...
%             'Marker', markers{k}, ...
%             'Color', colors(k,:), 'LineWidth', 1.5, ...
%             'MarkerSize', 6);
%     end
% end
% 
% title('Performance under varying Noise Levels (Single Surface)');
% xlabel('Noise Level (Standard Deviation)'); % Trục hoành bây giờ là mức nhiễu
% ylabel('RMSE (rad)');
% legend(algo_names, 'Location', 'northwest');
% set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
% xlim([min(noise_level_range), max(noise_level_range)]);



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
