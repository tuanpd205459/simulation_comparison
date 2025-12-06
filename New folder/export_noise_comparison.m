clc, clear, close all;
%%
% load("cchuan_so_sanh_cac_be_mat_nhieu.mat");
load("cchuan_ss_cac_thuat_toan_ban_cu.mat");
%% ===================== Configuration =====================
algoNames = {'Goldstein Branch-Cut', 'Quality-Guided', ...
             'Weighted Least Squares', 'Proposed Method'};
nAlgo     = numel(algoNames);

% Markers & color scheme (tối ưu hiển thị khoa học)
markers     = {'s','*','d','o'};
lineWidths  = 1.8;
markerSize  = 6;

% --- Màu (đảm bảo độ tương phản, phù hợp nền trắng) ---
colors = [
    0.75 0.00 0.00 ;   % 1. Goldstein - strong red
    0.96 0.47 0.15;   % 2. Quality-Guided - Bright Orange (high visibility)
    0.18 0.55 0.34;   % 3. WLS - Deep Green (natural, strong contrast)
    0.00 0.37 0.70    % 4. Proposed Method - Deep Science Blue (main highlight)
];

% ===================== Plot — RMSE Comparison =====================
figure('Color','w'); hold on;

for k = 1:nAlgo
    p(k) = plot(x_noise_axis, RMSE(:,k), ...
        'Color', colors(k,:), ...
        'LineWidth', lineWidths, ...
        'Marker', markers{k}, ...
        'MarkerSize', markerSize, ...
        'MarkerFaceColor', colors(k,:) * 0.9);                       % tạo contrast
end

% ===================== Labels & Legend =====================
% hTitle  = title('RMSE Comparison of Phase Unwrapping Algorithms');
hXLabel = xlabel('Noise Level (%)');
hYLabel = ylabel('RMSE (rad)');

hLegend = legend(p, algoNames, ...
    'Location','best', ...
    'FontSize',9, ...
    'Box','off');                   % tinh tế hơn

% ===================== Axes — Figure Formatting =====================
set(gca,'FontName','Helvetica', ...
        'FontSize',10, ...
        'LineWidth',1.2, ...
        'Box','off', ...
        'TickDir','out', ...
        'TickLength',[.015 .015], ...
        'XMinorTick','off','YMinorTick','off', ...
        'YGrid','on');

set(gcf, 'Units','Inches', 'Position', [1 1 6 4]);       % size figure (inch) 6x4

xlim([min(x_noise_axis)*1.05, max(x_noise_axis)*1.05]);
ylim([min(RMSE(:))- 0.025, max(RMSE(:)) * 1.05]); % Thiết lập giới hạn trục Y theo dữ liệu



% Đậm tiêu đề như MATLAB ví dụ
% set(hTitle,'FontSize',12);
set([hXLabel, hYLabel],'FontSize',11);

saveFolder = fullfile(pwd, 'so_sanh_cac_tt_nhieu_muc_nhieu');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

fileName = ['Fig_Comparison_noise' timestamp];   % đổi ten anh
fullPath = fullfile(saveFolder, fileName);
export_fig([fullPath '.png'], '-png', '-r600');       % PNG 600 dpi
export_fig([fullPath '.eps'], '-eps', '-opengl');   % EPS vector

% ===================== Export Figure =====================


% Tạo folder lưu kết quả
saveFolder = fullfile(pwd, 'so_sanh_cac_tt_nhieu_muc_nhieu');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Tên file có timestamp tự động
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
fileName = ['Fig_Comparison_noise' timestamp];   % đổi tên nếu muốn
fullPath = fullfile(saveFolder, fileName);

% Xuất PNG (600 dpi chuẩn journal)
print(gcf, [fullPath '.png'], '-dpng', '-r600');

% Xuất EPS vector (cho LaTeX/báo hội nghị)
print(gcf, [fullPath '.eps'], '-depsc', '-opengl');

disp("Ảnh đã lưu vào thư mục: " + saveFolder);
