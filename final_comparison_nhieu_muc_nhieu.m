clc, clear, close all;

load('cchuan_so_sanh_cac_be_mat_nhieu.mat');

%% --- XỬ LÝ DỮ LIỆU: ĐẨY GIÁ TRỊ VỀ DƯƠNG (OFFSET TỪ 0) ---
% Hàm nội bộ để đẩy min về 0
shift_to_zero = @(x) x - min(x(:));

% Áp dụng cho các biến (trừ Quality map vì thường nó đã là 0-1)
gt_display       = shift_to_zero(object_with_noise_cut);
wrapped_display  = shift_to_zero(wrapped_phase_cut);  
goldstein_display= shift_to_zero(phi_goldstein_cut);  
wls_display      = shift_to_zero(phi_wls_cut);          
proposed_display = shift_to_zero(phi_proposed_cut);
quality_display = shift_to_zero(phi_qual_cut);

%% --- FIGURE 1: HIỂN THỊ 6 HÌNH ẢNH (2 HÀNG x 3 CỘT) ---
figure('Name', 'Visual Comparison (Positive Offset)', 'Color', 'w', 'Position', [50, 50, 1200, 600]);
colormap turbo;

% Hàng 1
subplot(2, 3, 1);
imagesc(gt_display); axis image; axis off; colorbar;
title('(a)');

subplot(2, 3, 2);
imagesc(wrapped_phase_cut); axis image; axis off; colorbar;
title('(b)');

subplot(2, 3, 3);
imagesc(goldstein_display); axis image; axis off; colorbar;
title('(c)');

% Hàng 2
subplot(2, 3, 4);
imagesc(quality_display); axis image; axis off; colorbar;
title('(d)');

subplot(2, 3, 5);
imagesc(wls_display); axis image; axis off; colorbar;
title('(e)');

subplot(2, 3, 6);
imagesc(proposed_display); axis image; axis off; colorbar;
title('(f)');

% sgtitle('Phase Maps Comparison');
saveFolder = fullfile(pwd, 'so_sanh_cac_tt_nhieu_muc_nhieu');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

fileName = ['Fig_Comparison_ban2' timestamp];   % đổi ten anh
fullPath = fullfile(saveFolder, fileName);
export_fig([fullPath '.png'], '-png', '-r600');       % PNG 600 dpi
export_fig([fullPath '.eps'], '-eps', '-opengl');   % EPS vector
%% --- FIGURE 2: PHÂN TÍCH CHI TIẾT (MẶT CẮT & SAI SỐ) ---
figure('Name', 'Analysis: Cross-section & Error', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Lấy dòng giữa để cắt
mid_row = round(size(gt_display, 1) / 2);

% --- Subplot trên: MẶT CẮT NGANG (Cross-section) ---
subplot(1, 2, 1);
% Vẽ Ground Truth làm nền
plot(gt_display(mid_row, :), 'k', 'LineWidth', 2, 'DisplayName', 'GroundTruth'); hold on;
% Vẽ các thuật toán khác
plot(goldstein_display(mid_row, :), 'b--', 'LineWidth', 1, 'DisplayName', 'Goldstein');
plot(wls_display(mid_row, :),       'g:',  'LineWidth', 1.5, 'DisplayName', 'WLS');
plot(quality_display(mid_row, :),       'r:',  'LineWidth', 1.5, 'DisplayName', 'WLS');

plot(proposed_display(mid_row, :),  'c-.', 'LineWidth', 1.5, 'DisplayName', 'Proposed');
hold off;

% grid on; 
axis tight;
title(['Cross-section Profile at Row ', num2str(mid_row)]);
legend('Location', 'best');
ylabel('Phase Value');

% --- Subplot dưới: SAI SỐ (Error Plot) ---
subplot(1, 2, 2);
% Tên các thuật toán để hiển thị trong chú giải
ALGO_NAMES = {'Goldstein Branch-Cut', 'Quality-Guided', 'Weighted Least Squares', 'Proposed Method'};
N_ALGO = length(ALGO_NAMES);

% Cấu hình các đặc tính của đường vẽ (Marker, LineStyle, Color)
MARKERS = {'s', '*', 'd', 'o'}; % 's': Square, '*': Star, 'd': Diamond, 'o': Circle
LINE_STYLES = {'-', '-', '-', '-'}; % Nét liền cho tất cả
BASE_COLORS = lines(N_ALGO); % Mảng màu cơ bản

% Định nghĩa màu nổi bật cho thuật toán Proposed (thường là màu Cyan)
PROPOSED_COLOR = [0, 0.75, 0.75]; % Màu xanh ngọc (Cyan)

%% 2. Vẽ Biểu đồ So sánh RMSE theo Mức Nhiễu
% Tạo cửa sổ đồ thị
% figure('Name', 'RMSE vs Noise Level', 'Color', 'w', ...
%        'Position', [100, 100, 1000, 600]);

% --- Line Plot ---
hold on; 
% grid on; 
box on; % Thêm khung bao quanh đồ thị

% Vòng lặp vẽ dữ liệu cho từng thuật toán
for k = 1:N_ALGO
    % Xác định màu và cấu hình
    if strcmp(ALGO_NAMES{k}, 'Proposed Method')
        % Cấu hình nổi bật cho Proposed
        current_color = PROPOSED_COLOR;
        % marker_face_color = PROPOSED_COLOR;
        line_width = 1.5; % Nét dày hơn
    else
        % Cấu hình tiêu chuẩn cho các thuật toán còn lại
        current_color = BASE_COLORS(k,:);
        marker_face_color = 'none'; % Không tô màu mặt marker
        line_width = 1.5;
    end

    % Lệnh vẽ
    x_noise_axis = noise_level_range*100;
    plot(x_noise_axis, RMSE(:, k), ...
         'DisplayName', ALGO_NAMES{k}, ...
         'LineStyle', LINE_STYLES{k}, ...
         'Marker', MARKERS{k}, ...
         'Color', current_color, ...
         'LineWidth', line_width, ...
         'MarkerFaceColor', marker_face_color, ...
         'MarkerSize', 6); % Marker lớn hơn để dễ nhìn
end

% --- Định dạng Đồ thị ---
% Tiêu đề trục
xlabel('Noise Level (%)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('RMSE (rad)', 'FontSize', 12, 'FontName', 'Times New Roman');

% Chú giải
legend('Location', 'northwest', 'FontSize', 12, 'FontName', 'Times New Roman');

% Giới hạn trục x và y (tùy chọn)
xlim([min(x_noise_axis)*1.05, max(x_noise_axis)*1.05]);
ylim([min(RMSE(:))- 0.025, max(RMSE(:)) * 1.05]); % Thiết lập giới hạn trục Y theo dữ liệu

% Định dạng trục đồ thị
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', ...
         'XMinorGrid', 'off', 'YMinorGrid', 'off');
hold off;

saveFolder = fullfile(pwd, 'so_sanh_cac_tt_nhieu_muc_nhieu');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

fileName = ['Fig_Comparison_ban2' timestamp];   % đổi ten anh
fullPath = fullfile(saveFolder, fileName);
export_fig([fullPath '.png'], '-png', '-r600');       % PNG 600 dpi
export_fig([fullPath '.eps'], '-eps', '-opengl');   % EPS vector
