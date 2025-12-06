clc, clear, close all;
load("be_mat_khong_nhieu.mat");
%% Ve be mat 2D - turbo

addpath("D:\tuan\analysis\analysis_fringe\export_fig");
% Hàm nội bộ để đẩy min về 0
shift_to_zero = @(x) x - min(x(:));
groundtruth_cut = shift_to_zero(groundtruth_cut);
object_with_noise_cut       = shift_to_zero(object_with_noise_cut);
% wrapped_phase_cut  = shift_to_zero(wrapped_phase_cut);
phi_goldstein_cut= shift_to_zero(phi_goldstein_cut);
phi_wls_cut      = shift_to_zero(phi_wls_cut);
phi_proposed_cut = shift_to_zero(phi_proposed_cut);
phi_qual_cut = shift_to_zero(phi_qual_cut);
% --- Các thuật toán cần so sánh theo thu tu chuan---
dataList = { groundtruth_cut,...
    wrapped_phase_cut,...
    phi_goldstein_cut, ...
    phi_qual_cut, ...
    phi_wls_cut, ...
    phi_proposed_cut };
phase_names = {'Goldstein', ...
    'Quality','2D-WLS', ...
    'Proposed'};
nPhase = length(dataList);


figure;
surf(phi_wls_cut,"EdgeColor","none");

% 1.5 SETUP SPATIAL COORDINATES (MM)
[rows, cols] = size(dataList{1,1});
% px_size = 3.45e-3; % 3.45 µm = 0.00345 mm
% x_vec = (0 : cols-1) * px_size;
% y_vec = (0 : rows-1) * px_size;
% giu nguyen truc toa do la pixel
x_vec = (0 : cols-1) ;
y_vec = (0 : rows-1) ;
%
all_pixels = []; 
for i = 1:size(dataList, 1)
    d = dataList{i,1};
    all_pixels = [all_pixels; d(:)]; 
end
robust_min = prctile(all_pixels, 0.1); 
robust_max = prctile(all_pixels, 99.8); 
z_lims = [robust_min, robust_max];
clear all_pixels;
%% 3. FIGURE SETTINGS
figWidth  = 17.5;
figHeight = 10;
fontSize  = 10;
fontName  = 'Times New Roman';

fig = figure('Units', 'centimeters', ...
             'Position', [2, 2, figWidth, figHeight], ...
             'Color', 'w', ...
             'Name', 'Fig_Comparison_5_simulation_2D_MM_turbo', ...
             'NumberTitle', 'off');

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 4. DRAW 5 SUBFIGURES (2D)
num_imgs = length(dataList);
cols_fig = 3;     % số cột của layout

labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
axs = gobjects(1, length(dataList));

for i = 1: length(dataList)
    axs(i) = nexttile;

    data = dataList{i,1};

    imagesc(x_vec, y_vec, data);
    axis image;
    clim(z_lims);
    colormap(gca, turbo);

    % Thêm nhãn (a), (b), ...
    title(labels{i}, 'FontWeight','bold', 'FontSize', fontSize+1, ...
        'FontName','Times New Roman', 'Interpreter','latex');

    xlabel('x (pixel)', 'Interpreter', 'latex');
    ylabel('y (pixel)', 'Interpreter', 'latex');

    set(gca, 'FontName', fontName, 'FontSize', fontSize, ...
        'LineWidth', 1, 'TickLabelInterpreter', 'latex');
    box on;
end
cb = colorbar;
cb.Layout.Tile = 'east'; 
cb.Limits = z_lims;

cb.TickLabelInterpreter = 'latex';
cb.FontSize = fontSize;

cb.Label.String = 'Phase (rad)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fontSize + 1;

% saveFolder = fullfile(pwd, 'ExportedFigures_experiments');
% if ~exist(saveFolder, 'dir')
%     mkdir(saveFolder);
% end
% timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% 
% fileName = ['Fig_Comparison_5_Real_2D_MM_turbo' timestamp];   % đổi ten anh
% fullPath = fullfile(saveFolder, fileName);
% export_fig([fullPath '.png'], '-png', '-r600');       % PNG 600 dpi
% export_fig([fullPath '.eps'], '-eps', '-opengl');   % EPS vector

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
