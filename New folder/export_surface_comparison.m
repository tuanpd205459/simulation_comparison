clc, clear, close all;
load("cchuan_ss_cac_be_mat.mat");
addpath("D:\tuan\analysis\analysis_fringe\export_fig");

%% ===================== Configuration - line + box plot
%% ===================== Configuration =====================
algoNames = {'Goldstein Branch-Cut', 'Quality-Guided', ...
             'Weighted Least Squares', 'Proposed Method'};
nAlgo     = numel(algoNames);

markers     = {'s','*','d','o'};
markerSize  = 6;
lineWidth   = 1.8;

% ================= Color Palette (Scientific Professional) =================
colors = [
    0.75 0.00 0.00 ;   % Goldstein - strong red
    0.96 0.47 0.15;   % Quality - Bright Orange
    0.18 0.55 0.34;   % WLS - Deep Green
    0.00 0.37 0.70    % Proposed - Deep Science Blue (highlight)
    ];

%% ===================== Create Figure =====================
figure('Color','w','Position',[200 100 900 380]);

%% ===================== (a) Line Plot =====================
subplot(1,2,1); hold on;
for k = 1:nAlgo
    p(k) = plot(1:nRuns, RMSE(:,k), ...
        'Color', colors(k,:), ...
        'LineWidth', lineWidth, ...
        'Marker', markers{k}, ...
        'MarkerSize', markerSize, ...
        'MarkerFaceColor', colors(k,:) ...
        );      % tạo contrast khi in grayscale
end

title('(a)');
xlabel('Surface Index','FontSize',11);
ylabel('RMSE (rad)','FontSize',11);

legend(p, algoNames,'FontSize',9,'Location','best','Box','off');
xlim([1 nRuns]);

% Axes style
set(gca,'FontName','Helvetica','FontSize',10,'LineWidth',1.2,...
        'Box','off','TickDir','out','XMinorTick','off','YMinorTick','off',...
        'TickLength',[.015 .015],'YGrid','on');

%% ===================== (b) Boxplot =====================
subplot(1,2,2);
boxplot(RMSE,'Labels',algoNames,'Symbol','');

title('(b)');
ylabel('RMSE (rad)','FontSize',11);

set(gca,'FontName','Helvetica','FontSize',10,'LineWidth',1.2,...
        'Box','off','TickDir','out','YGrid','on','TickLength',[.015 .015]);

% --------- Fill Color for Boxplot ---------
h = findobj(gca,'Tag','Box');
for k = 1:nAlgo
    patch(get(h(nAlgo-k+1),'XData'), get(h(nAlgo-k+1),'YData'), colors(k,:), ...
          'FaceAlpha',0.35,'EdgeColor','k','LineWidth',1.2);
end
% Fix legend/axes bị cắt khi export
set(gca,'LooseInset',max(get(gca,'TightInset')*1.3,0.05));
%%
% saveFolder = fullfile(pwd, 'so_sanh_nhieu_be_mat');
% if ~exist(saveFolder, 'dir')
%     mkdir(saveFolder);
% end
% timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% 
% fileName = ['Fig_Comparison_nhieu_mat' timestamp];   % đổi ten anh
% fullPath = fullfile(saveFolder, fileName);
% export_fig([fullPath '.png'], '-png', '-r600');       % PNG 600 dpi
% export_fig([fullPath '.eps'], '-eps', '-opengl');   % EPS vector
%% ===================== Export Figure =====================

% Tạo folder lưu kết quả
saveFolder = fullfile(pwd, 'so_sanh_nhieu_be_mat');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Tên file có timestamp tự động
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
fileName = ['Fig_Comparison_nhieu_mat_' timestamp];   % đổi tên nếu muốn
fullPath = fullfile(saveFolder, fileName);

% Xuất PNG (600 dpi chuẩn journal)
print(gcf, [fullPath '.png'], '-dpng', '-r600');

% Xuất EPS vector (cho LaTeX/báo hội nghị)
print(gcf, [fullPath '.eps'], '-depsc', '-opengl');

disp("Ảnh đã lưu vào thư mục: " + saveFolder);


%%

%% ===================== Configuration =====================
algoNames = {'Goldstein Branch-Cut', 'Quality-Guided', ...
             'Weighted Least Squares', 'Proposed Method'};
nAlgo     = numel(algoNames);

% Markers & color scheme (tối ưu hiển thị khoa học)
markers     = {'s','*','d','o'};
lineWidths  = 1.5;
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
    p(k) = plot(1:nRuns, RMSE(:,k), ...
        'Color', colors(k,:), ...
        'LineWidth', lineWidths, ...
        'Marker', markers{k}, ...
        'MarkerSize', markerSize, ...
        'MarkerFaceColor', colors(k,:) * 0.9);                       % tạo contrast
end

% ===================== Labels & Legend =====================
% hTitle  = title('RMSE Comparison of Phase Unwrapping Algorithms');
hXLabel = xlabel('Surface Index');
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

xlim([1 nRuns]);
% ylim auto tùy dữ liệu, nếu muốn:
ylim([0 max(RMSE(:))*1.1])

% Đậm tiêu đề như MATLAB ví dụ
% set(hTitle,'FontSize',12,'FontWeight','bold');
set([hXLabel, hYLabel],'FontSize',11);

saveFolder = fullfile(pwd, 'so_sanh_nhieu_be_mat');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

fileName = ['Fig_Comparison_nhieu_mat' timestamp];   % đổi ten anh
fullPath = fullfile(saveFolder, fileName);
export_fig([fullPath '.png'], '-png', '-r600');       % PNG 600 dpi
export_fig([fullPath '.eps'], '-eps', '-opengl');   % EPS vector

%% ===================== Export Figure =====================


% Tạo folder lưu kết quả
saveFolder = fullfile(pwd, 'so_sanh_nhieu_be_mat');
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Tên file có timestamp tự động
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
fileName = ['Fig_Comparison_nhieu_mat_' timestamp];   % đổi tên nếu muốn
fullPath = fullfile(saveFolder, fileName);

% Xuất PNG (600 dpi chuẩn journal)
print(gcf, [fullPath '.png'], '-dpng', '-r600');

% Xuất EPS vector (cho LaTeX/báo hội nghị)
print(gcf, [fullPath '.eps'], '-depsc', '-opengl');

disp("Ảnh đã lưu vào thư mục: " + saveFolder);
