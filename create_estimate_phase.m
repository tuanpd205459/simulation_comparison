function [estimate_phase, wrapped_phase] = create_estimate_phase(hologram, groundtruth)
object_phase = groundtruth;
SENSITIVITY_COEF = 0.6;   % Hệ số nhạy adaptive threshold
NEIGHBORHOOD_SIZE = 51;    % Kích thước vùng lân cận (phải là số lẻ)
MIN_BRANCH_LENGTH = 8;    % Ngưỡng độ dài nhánh thừa cần xóa (thay cho distThresh cũ)
GAUSS_SIGMA = 1;           % Độ làm mượt ảnh
distThresh = 6;

hologram = imgaussfilt(hologram, GAUSS_SIGMA);
hologram = adapthisteq(hologram);

T = adaptthresh(hologram, SENSITIVITY_COEF, ...
    'NeighborhoodSize', [NEIGHBORHOOD_SIZE NEIGHBORHOOD_SIZE], ...
    'Statistic', 'median');

BW = imbinarize(hologram, T);

BW_skel_clean = bwskel(BW, 'MinBranchLength', MIN_BRANCH_LENGTH);

BW_skel_clean = bwmorph(BW_skel_clean, 'clean');

% cắt cách điểm endpoint và branchpoints gần nhau
BW = BW_skel_clean;
MIN_BRANCH_LENGTH = 10;
BW = bwskel(BW, 'MinBranchLength', MIN_BRANCH_LENGTH);

%%
BW = bwfill(BW,'holes');
BW = bwskel(BW, 'MinBranchLength', MIN_BRANCH_LENGTH);

%% ngawts ket noi H-brigde
B = bwmorph(BW, 'branchpoints');
R = 4;
se_disk = strel('disk', R);

B = imdilate(B, se_disk);
E = bwmorph(BW, 'endpoints');
BW = BW & ~B;
min_len = 30; % Ngưỡng độ dài (pixel)
BW = bwareaopen(BW, min_len);
BW = BW | B;
BW = bwmorph(BW, 'thin', Inf);

%% --- Tìm endpoint ---

BW = bwmorph(BW,"bridge",Inf);
BW = bwmorph(BW,"diag", Inf);
BW = bwmorph(BW,"skeleton", Inf);
BW = bwmorph(BW,'spur',1);

%% noois diem
n = 8;
for count =1:n
    % Timf endPoint
    % Kernel để đếm số hàng xóm (8-neighbors)
    endPoints = bwmorph(BW, 'endpoints');
    vectors = fitEndpointVectors(BW, endPoints, 20);
    %
    if count == 1 %nối vân dài + góc lệch nhỏ + khoảng cách nhỏ
        minCompSize = 12;   % chỉ nối nếu component đủ dài
        maxDist     = 6;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(15);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count ==2 % Nối vân dài + góc lệch lớn + khoảng cách lớn
        minCompSize = 12;   % chỉ nối nếu component đủ dài
        maxDist     = 12;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(15);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count == 3 % Nối vân dài + góc lệch nhỏ + khoảng cách lớn hơn
        minCompSize = 12;   % chỉ nối nếu component đủ dài
        maxDist     = 25;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(30);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count == 4 % Nối vân ngắn + góc lệch lớn
        minCompSize = 5;   % chỉ nối nếu component đủ dài
        maxDist     = 50;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(40);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count == 5 % Nối vân ngắn + góc lệch lớn
        minCompSize = 20;   % chỉ nối nếu component đủ dài
        maxDist     = 50;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(15);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count == 6 % Nối vân ngắn + góc lệch lớn
        minCompSize = 20;   % chỉ nối nếu component đủ dài
        maxDist     = 50;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(30);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count == 7 % Nối vân ngắn + góc lệch lớn
        minCompSize = 20;   % chỉ nối nếu component đủ dài
        maxDist     = 50;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(30);  % = 0.866 ~ hướng lệch <= 30°
    end
    if count == 8 % Nối vân ngắn + góc lệch lớn
        minCompSize = 20;   % chỉ nối nếu component đủ dài
        maxDist     = 60;   % khoảng cách tối đa giữa 2 endpoint
        vecAlignThr = cosd(80);  % = 0.866 ~ hướng lệch <= 30°
    end
    % BW = bwmorph(BW,"spur",1);
    CC = bwconncomp(BW, 8);
    [BW, linesConnected] = connectEndpoints(BW, vectors, CC, minCompSize, maxDist, vecAlignThr);
end

%% Nối ở biên
% --- Tìm endpoint ---
endPoints = bwmorph(BW,'endpoints');
vectors = fitEndpointVectors(BW, endPoints, 12);
margin = 30;
extendLength = 30;
% --- Tham số nối ---
BW = extendLineNearBorder(BW, vectors, extendLength, margin);

endPoints_mask = bwmorph(BW, 'endpoints');
vectors_data = fitEndpointVectors(BW, endPoints_mask, 12);

if ~isempty(vectors_data)
    % Tách toạ độ và vector hướng từ kết quả tính được
    endpoints_ordinate = vectors_data(:, 1:2); % Cột 1,2 là x,y
    vectors_dir      = vectors_data(:, 3:4);   % Cột 3,4 là vx,vy

    % 2. Cấu hình tham số nối
    maxLen_final = 50; % Độ dài tìm kiếm (tăng lên nếu vân đứt xa)

    % 3. Gọi hàm nối giao cắt (Hàm mới sửa dùng CCW)
    [BW_connected, list_con] = connect_intersecting_ridges(BW, endpoints_ordinate, vectors_dir, maxLen_final);

    % Vẽ các đoạn vừa nối thêm để dễ nhìn
    if ~isempty(list_con)
        for i = 1:size(list_con, 1)
            % list_con chứa [P1_x P1_y P2_x P2_y] hoặc [P1; P2] tuỳ cách hàm trả về
            % Ở code trước trả về [P1, P2] tức là 1 hàng có 4 phần tử
            plot([list_con(i,1) list_con(i,3)], [list_con(i,2) list_con(i,4)], 'g-', 'LineWidth', 2);
        end
    end
else
    BW_connected = BW;
end
MIN_BRANCH_LENGTH = 40;
BW_connected = bwskel(BW_connected, 'MinBranchLength', MIN_BRANCH_LENGTH);

%%
MIN_BRANCH_LENGTH = 10;
BW_connected = bwskel(BW_connected, 'MinBranchLength', MIN_BRANCH_LENGTH);
offset = 15;
BW_connected = BW_connected(offset : end - offset+ 1, ...
    offset : end - offset + 1);
branchpoints_check = bwmorph(BW_connected, "branchpoints");
% kiem tra neu con branchpoint thi dung
if any(branchpoints_check(:))
    % Thực hiện lệnh dừng hoặc báo lỗi tại đây
    disp('Phát hiện điểm rẽ nhánh! Đang dừng...');
    return; % Hoặc break
end
%%
wrapped_phase = wrapToPi(object_phase) ;

%% tái tạo pha estimate
[~, labels, img] = assign_fringe_order(BW_connected, true);

% Tái tạo bề mặt từ vân
[phi_est, ~] = reconSurface_linearPushed(img, labels, 632.8e-9, 'None', false);
phi_est = phi_est(5:end-5, 5:end-5);
phi_est = phi_est - min(phi_est(:));

[N,M] = size(groundtruth);

fx = 20 / N; % Tần số sóng mang
fy = -30 / M;
[X, Y] = meshgrid(1:N, 1:M);
plane_phase = 2*pi*(fx*X + fy*Y);

plane_phase = plane_phase - min(plane_phase(:));
[phi_est, plane_phase,...
    wrapped_phase] = crop_multiple_to_smallest(phi_est, plane_phase,...
   wrapped_phase);
phi_est = phi_est - plane_phase - (max(phi_est(:)- max(plane_phase(:))))/2;

[ phi_est,wrapped_phase]...
    = crop_multiple_to_smallest(phi_est, wrapped_phase);
[finalUnwrappedPhase, kMap] = unwrapUsingEstimate(phi_est, wrapped_phase);

offset = 7;
finalUnwrappedPhase = finalUnwrappedPhase(offset+1:end-offset, offset+1:end-offset);
estimate_phase = finalUnwrappedPhase;
end


%%
function [BW_new, linesConnected] = connectEndpoints(BW, vectors, CC, minCompSize, maxDist, vecAlignThr)
% BW           : skeleton binary
% vectors      : [cx cy vx vy] từ hàm computeEndpointVectors
% CC           : bwconncomp(BW,8)
% minCompSize  : kích thước tối thiểu của vân
% maxDist      : khoảng cách tối đa cho phép nối
% vecAlignThr  : ngưỡng cos(angle) hướng (ví dụ 0.7 ~ >45°)

BW_new = BW; % copy để cập nhật nối
linesConnected = {}; % cell lưu danh sách các đoạn đã nối

for i = 1:size(vectors,1)-1
    cx1 = vectors(i,1); cy1 = vectors(i,2);
    v1  = [vectors(i,3), vectors(i,4)];

    % kiểm tra component của endpoint i
    comp_id1 = findComponent(CC, [cy1,cx1]);
    if comp_id1==0 || numel(CC.PixelIdxList{comp_id1}) < minCompSize
        continue;
    end

    for j = i+1:size(vectors,1)
        cx2 = vectors(j,1); cy2 = vectors(j,2);
        v2  = [vectors(j,3), vectors(j,4)];

        % kiểm tra component j
        comp_id2 = findComponent(CC, [cy2,cx2]);
        if comp_id2==0 || numel(CC.PixelIdxList{comp_id2}) < minCompSize
            continue;
        end

        % --- khoảng cách Euclidean giữa 2 endpoint ---
        d = hypot(cx1-cx2, cy1-cy2);
        if d > maxDist, continue; end

        % --- kiểm tra hướng vector (cùng hướng nối) ---
        dir12 = [cx2-cx1, cy2-cy1];
        dir12 = dir12 / (norm(dir12)+eps);

        cond1 = dot(v1, dir12) > vecAlignThr;    % v1 hướng về P2
        cond2 = dot(v2, -dir12) > vecAlignThr;   % v2 hướng về P1

        if ~(cond1 && cond2), continue; end

        % --- kiểm tra thêm khoảng cách vuông góc ---
        % đường thẳng qua P2 với vector v2
        a = -v2(2);
        b =  v2(1);
        c =  v2(2)*cx2 - v2(1)*cy2;
        d_perp = abs(a*cx1 + b*cy1 + c) / sqrt(a^2 + b^2);

        if d_perp > 5, continue; end

        % --- nối 2 endpoint ---
        [BW_new, linePixels] = drawLine(BW_new, cx1, cy1, cx2, cy2);
        linesConnected{end+1} = linePixels; %#ok<AGROW>
    end
end
end
function vectors = fitEndpointVectors(BW, endPoints, Nfit)
% fitEndpointVectors - Tính vector hướng tại endpoint của skeleton
%
% Cú pháp:
%   vectors = fitEndpointVectors(BW, endPoints, Nfit)
%
% Input:
%   BW        - ảnh nhị phân skeleton
%   endPoints - ảnh nhị phân endpoint (1 tại endpoint)
%   Nfit      - số pixel dùng để fit PCA (ví dụ: 30)
%
% Output:
%   vectors - ma trận [N x 4], mỗi hàng:
%             [cx cy vx vy]
%             (cx, cy) = tọa độ endpoint
%             (vx, vy) = vector đơn vị hướng ra ngoài

[y_idx, x_idx] = find(endPoints);  % tọa độ endpoints
CC = bwconncomp(BW, 8);           % tìm các component
vectors = [];

for k = 1:length(x_idx)
    cx = x_idx(k);
    cy = y_idx(k);

    % Kiểm tra endpoint thuộc component nào
    comp_id = 0;
    for c = 1:CC.NumObjects
        if ismember(sub2ind(size(BW), cy, cx), CC.PixelIdxList{c})
            comp_id = c;
            break;
        end
    end

    if comp_id == 0, continue; end  % endpoint không thuộc component nào

    % Lấy tọa độ tất cả pixel trong component
    [yy, xx] = ind2sub(size(BW), CC.PixelIdxList{comp_id});

    % Tính khoảng cách từ endpoint
    dist2 = (xx - cx).^2 + (yy - cy).^2;
    [~, idx] = sort(dist2);
    idxN = idx(1:min(Nfit, numel(idx)));

    X = xx(idxN);
    Y = yy(idxN);

    if numel(X) > 1
        % --- Fit hướng bằng PCA ---
        Xc = X - mean(X);
        Yc = Y - mean(Y);
        D = [Xc(:) Yc(:)];
        [~,~,V] = svd(D,'econ');
        v = V(:,1);  % vector chính (cột đầu tiên)
        v = v / norm(v);

        % --- Xác định hướng "ra ngoài" ---
        centroid = [mean(X); mean(Y)];
        c = centroid - [cx; cy];  % vector từ endpoint vào trong component
        if dot(v, c) > 0
            v = -v; % đảo dấu để hướng ra ngoài
        end
    else
        v = [0;0];
    end

    vectors = [vectors; cx cy v(1) v(2)];
end
end
function comp_id = findComponent(CC, p)
% p = [row, col]
comp_id = 0;
idx = sub2ind(CC.ImageSize, p(1), p(2));
for c = 1:CC.NumObjects
    if ismember(idx, CC.PixelIdxList{c})
        comp_id = c;
        return;
    end
end
end

function [BW, linePix] = drawLine(BW, x1, y1, x2, y2)
% Vẽ line nối từ (x1,y1) đến (x2,y2) bằng thuật toán Bresenham
[h, w] = size(BW);
[lineX, lineY] = bresenham(x1, y1, x2, y2);

linePix = [lineX(:), lineY(:)];

for k = 1:length(lineX)
    cx = lineX(k);
    cy = lineY(k);
    if cx >= 1 && cx <= w && cy >= 1 && cy <= h
        BW(cy, cx) = 1;
    end
end

end

function [x, y] = bresenham(x1, y1, x2, y2)

% Thuật toán Bresenham
x1 = round(x1); y1 = round(y1);
x2 = round(x2); y2 = round(y2);

dx = abs(x2 - x1);
dy = abs(y2 - y1);

sx = sign(x2 - x1);
sy = sign(y2 - y1);

err = dx - dy;

x = []; y = [];
while true
    x(end+1) = x1;
    y(end+1) = y1;
    if x1 == x2 && y1 == y2
        break;
    end
    e2 = 2 * err;
    if e2 > -dy
        err = err - dy;
        x1 = x1 + sx;
    end
    if e2 < dx
        err = err + dx;
        y1 = y1 + sy;
    end
end
end
%%

%% ========== Hàm chính ==========
function [BW_out, connections] = extend_and_connect(BW, endpoints, vectors, maxLen, step, connectThresh)
% EXTEND_AND_CONNECT
%   Extend rays from endpoints along vectors, then connect endpoints
%   whose rays meet or get close.
%
% INPUT:
%   BW         - ảnh nhị phân ban đầu
%   endpoints  - Nx2 [x y]
%   vectors    - Nx2 hướng tương ứng
%   maxLen     - chiều dài mở rộng
%   step       - bước nội suy
%   connectThresh - ngưỡng để nối (pixel)
%
% OUTPUT:
%   BW_out     - ảnh sau khi nối
%   connections - cell M×2 chứa cặp điểm được nối

assert(size(endpoints,2)==2, 'endpoints must be Nx2 [x y]');
assert(size(vectors,1)==size(endpoints,1), 'vectors must match endpoints count');

[H,W] = size(BW);
N = size(endpoints,1);

rays = cell(N,1);

% ---- 1) TẠO CÁC RAY DỌC THEO VECTOR ----
for i = 1:N
    p0 = endpoints(i,:);        % 1×2
    v  = vectors(i,:);          % 1×2

    if all(v == 0)
        rays{i} = p0;
        continue;
    end

    vn = v / norm(v);           % chuẩn hóa 1×2
    ts = (0:step:maxLen)';      % L×1

    % Tạo L×2 bằng nhân ma trận: ts*(1×2) hợp lệ → (L×2)
    pts = ts * vn + p0;         % L×2

    % Giữ điểm trong ảnh
    valid = pts(:,1)>=1 & pts(:,1)<=W & pts(:,2)>=1 & pts(:,2)<=H;
    pts = round(pts(valid,:));

    % tránh trùng pixel
    pts = unique(pts,'rows','stable');

    rays{i} = pts;
end

% ---- 2) TÌM CẶP RAY NÀO GẦN NHAU ----
connections = {};
for i = 1:N-1
    Pi = double(rays{i});
    if isempty(Pi), continue; end

    for j = i+1 : N
        Pj = double(rays{j});
        if isempty(Pj), continue; end

        % Tính khoảng cách qua pdist2
        D = pdist2(Pi, Pj);
        [dmin, idx] = min(D(:));

        if dmin <= connectThresh
            [ia, jb] = ind2sub(size(D), idx);

            p_closest = Pi(ia,:);
            q_closest = Pj(jb,:);

            connections(end+1,1:2) = {p_closest, q_closest};
        end
    end
end

% ---- 3) VẼ ĐƯỜNG NỐI TRÊN ẢNH ----
BW_out = BW;
for k=1:size(connections,1)
    BW_out = drawLineBW(BW_out, connections{k,1}, connections{k,2});
end
end

function [BW_new, connections] = connect_intersecting_ridges(BW, endpoints, vectors, maxLen)
BW_new = BW;
connections = [];

if isempty(endpoints), return; end
num_pts = size(endpoints, 1);

% 1. Định danh vùng liên thông
[L, ~] = bwlabel(BW);
ind = sub2ind(size(BW), round(endpoints(:,2)), round(endpoints(:,1)));
pt_labels = L(ind);

% Chuẩn hóa vector
norms = sqrt(sum(vectors.^2, 2));
valid = norms > 0;
vectors(valid, :) = bsxfun(@rdivide, vectors(valid, :), norms(valid));
vectors(~valid, :) = 0;

scale_factor = maxLen * 1.5;

for i = 1:num_pts
    for j = i+1:num_pts

        % Nếu cùng thuộc 1 vân thì bỏ qua
        if pt_labels(i) == pt_labels(j) && pt_labels(i) > 0
            continue;
        end

        P1 = endpoints(i, :);
        P2 = endpoints(j, :);

        % Bỏ qua nếu 2 đầu mút quá xa nhau (Filter sơ bộ)
        if norm(P1 - P2) > maxLen * 2
            continue;
        end

        V1 = vectors(i, :);
        V2 = vectors(j, :);

        P1_end = P1 + V1 * scale_factor;
        P2_end = P2 + V2 * scale_factor;

        connected = false;

        % --- TRƯỜNG HỢP 1: Hai tia cắt nhau hình học (Giao thoa giữa hư không) ---
        is_cut = segments_intersect(P1, P1_end, P2, P2_end);
        if is_cut
            if norm(P1 - P2) < maxLen * 1.5
                [BW_new, ~] = drawLine(BW_new, P1(1), P1(2), P2(1), P2(2));
                connections = [connections; P1, P2];
                connected = true;
            end
        end

        % --- TRƯỜNG HỢP 2: Ray Casting (Bắn tia tìm vùng) ---
        % Logic mới: Điểm trúng (Hit Point) phải gần endpoint đích

        % Kiểm tra tia P1 -> bắn về phía P2 (nhưng trúng vùng chứa P2)
        if ~connected
            hit_point_1 = trace_ray_to_component(P1, P1_end, L, pt_labels(j));

            if ~isempty(hit_point_1)
                dist_gap = norm(P1 - hit_point_1);      % Độ dài đoạn nối
                dist_to_tip = norm(hit_point_1 - P2);   % Độ lệch so với endpoint đích

                % ĐIỀU KIỆN QUAN TRỌNG MỚI THÊM:
                % dist_to_tip < maxLen: Nghĩa là điểm va chạm không được nằm quá xa đầu mút P2
                % Tránh việc P1 đâm vào "bụng" của vân chứa P2
                if dist_gap < maxLen * 1.5 && dist_to_tip < maxLen
                    [BW_new, ~] = drawLine(BW_new, P1(1), P1(2), hit_point_1(1), hit_point_1(2));
                    connections = [connections; P1, hit_point_1];
                    connected = true;
                end
            end
        end

        % Kiểm tra tia P2 -> bắn về phía P1
        if ~connected
            hit_point_2 = trace_ray_to_component(P2, P2_end, L, pt_labels(i));

            if ~isempty(hit_point_2)
                dist_gap = norm(P2 - hit_point_2);
                dist_to_tip = norm(hit_point_2 - P1); % So với đầu mút P1

                if dist_gap < maxLen * 1.5 && dist_to_tip < maxLen
                    [BW_new, ~] = drawLine(BW_new, P2(1), P2(2), hit_point_2(1), hit_point_2(2));
                    connections = [connections; P2, hit_point_2];
                end
            end
        end

    end
end

if size(BW_new, 3) > 1
    BW_new = imbinarize(rgb2gray(BW_new));
end
end

% --- Giữ nguyên hàm phụ trợ ---
function hit_point = trace_ray_to_component(P_start, P_end, LabelMatrix, target_label)
hit_point = [];
if target_label == 0, return; end

dist = norm(P_end - P_start);
num_steps = ceil(dist);

x_vals = linspace(P_start(1), P_end(1), num_steps);
y_vals = linspace(P_start(2), P_end(2), num_steps);

[H, W] = size(LabelMatrix);
start_offset = 3;

for k = start_offset:num_steps
    cx = round(x_vals(k));
    cy = round(y_vals(k));

    if cx < 1 || cx > W || cy < 1 || cy > H, break; end

    pixel_label = LabelMatrix(cy, cx);

    if pixel_label == target_label
        hit_point = [cx, cy];
        return;
    elseif pixel_label > 0 && pixel_label ~= target_label
        return;
    end
end
end

% Hàm kiểm tra cắt nhau
function res = segments_intersect(p1, p2, p3, p4)
d1 = ccw(p3, p4, p1); d2 = ccw(p3, p4, p2);
d3 = ccw(p1, p2, p3); d4 = ccw(p1, p2, p4);
res = ((d1 * d2) < 0) && ((d3 * d4) < 0);
end
function val = ccw(a, b, c)
val = (b(1) - a(1)) * (c(2) - a(2)) - (b(2) - a(2)) * (c(1) - a(1));
end
% --- Hàm phụ trợ: Bắn tia tìm vùng liên thông đích ---

function BW_out = extendLineNearBorder(BW, vectors, extendLen, margin)
% extendLineNearBorder - Nối dài endpoint ra ngoài NẾU nó gần biên ảnh
%
% Input:
%   BW        - ảnh nhị phân
%   vectors   - [cx, cy, vx, vy] cho mỗi endpoint
%   extendLen - số pixel muốn nối dài thêm
%   margin    - ngưỡng khoảng cách từ biên (ví dụ 5)
%
% Output:
%   BW_out    - ảnh nhị phân sau khi vẽ đoạn thẳng nối dài

[H,W] = size(BW);
BW_out = BW;

for i = 1:size(vectors,1)
    cx = vectors(i,1);
    cy = vectors(i,2);
    vx = vectors(i,3);
    vy = vectors(i,4);

    % --- CHỈ vẽ nếu endpoint gần biên ---
    if ~(cx <= margin || cx >= W-margin || cy <= margin || cy >= H-margin)
        continue; % bỏ qua nếu không gần biên
    end

    % Tính điểm mới (C = B + extendLen*v)
    x3 = cx + extendLen*vx;
    y3 = cy + extendLen*vy;

    % Bresenham từ (cx,cy) đến (x3,y3)
    [xLine, yLine] = bresenham2(round(cx), round(cy), round(x3), round(y3));

    % Loại pixel ngoài biên
    mask = xLine>=1 & xLine<=W & yLine>=1 & yLine<=H;
    xLine = xLine(mask);
    yLine = yLine(mask);

    % Vẽ vào ảnh
    BW_out(sub2ind([H,W], yLine, xLine)) = 1;
end

end
function [x,y] = bresenham2(x1,y1,x2,y2)
x1=round(x1); y1=round(y1);
x2=round(x2); y2=round(y2);

dx=abs(x2-x1); dy=abs(y2-y1);
sx=sign(x2-x1); sy=sign(y2-y1);

x=x1; y=y1;
xx=[]; yy=[];

if dx > dy
    err = dx/2;
    while x ~= x2
        xx(end+1)=x; yy(end+1)=y;
        x = x + sx;
        err = err - dy;
        if err < 0
            y = y + sy;
            err = err + dx;
        end
    end
else
    err = dy/2;
    while y ~= y2
        xx(end+1)=x; yy(end+1)=y;
        y = y + sy;
        err = err - dx;
        if err < 0
            x = x + sx;
            err = err + dy;
        end
    end
end
xx(end+1)=x2; yy(end+1)=y2;
x=xx; y=yy;
end

function [fringe_order, fringe_labels, processed_image] = assign_fringe_order(input_image, display_result)
% ASSIGN_FRINGE_ORDER Gán bậc vân cho ảnh hologram đã được skeletonize
%
% Hàm này thực hiện gán nhãn bậc vân dựa trên vị trí tương đối so với tâm ảnh.
% Vân gần tâm nhất được gán bậc 0, các vân phía trên có bậc dương tăng dần,
% các vân phía dưới có bậc âm giảm dần.
%
% INPUT:
%   input_image    - Ảnh binary đã được skeletonize
%   display_result - (Optional) true/false để hiển thị kết quả (default: true)
%
% OUTPUT:
%   fringe_order     - Số lượng vân được phát hiện
%   fringe_labels    - Vector chứa nhãn bậc vân của từng vùng liên thông
%   processed_image  - Ảnh đã được cắt biên và xử lý
%
% EXAMPLE:
%   [order, labels, img] = assign_fringe_order(skeleton_image);
%   [order, labels, img] = assign_fringe_order(skeleton_image, false); % Không hiển thị

% --- Xử lý tham số đầu vào ---
if nargin < 1
    error('Thiếu tham số đầu vào: input_image');
end

if nargin < 2
    display_result = true; % Mặc định hiển thị kết quả
end

% --- Kiểm tra đầu vào ---
if isempty(input_image)
    error('Ảnh đầu vào không được để trống');
end

if ~islogical(input_image) && ~(isnumeric(input_image) && all(input_image(:) == 0 | input_image(:) == 1))
    error('Ảnh đầu vào phải là ảnh binary (logical hoặc 0/1)');
end

% Chuyển đổi sang logical nếu cần
if ~islogical(input_image)
    input_image = logical(input_image);
end

try
    % --- Bước 1: Cắt biên ảnh để tránh ảnh hưởng vùng biên ---
    offset = 0;
    [orig_H, orig_W] = size(input_image);

    % Kiểm tra kích thước ảnh
    if orig_H <= 2*offset || orig_W <= 2*offset
        warning('Ảnh quá nhỏ để cắt biên. Sử dụng ảnh gốc.');
        bw_crop = input_image;
        offset = 0;
    else
        bw_crop = input_image(offset+1:end-offset, offset+1:end-offset);
    end

    [H, W] = size(bw_crop);

    % --- Bước 2: Tìm các vùng liên thông (vân) ---

    cc = bwconncomp(bw_crop);

    if cc.NumObjects == 0
        warning('Không tìm thấy vân nào trong ảnh');
        fringe_order = 0;
        fringe_labels = [];
        processed_image = bw_crop;
        return;
    end

    labeled_matrix = labelmatrix(cc);
    stats = regionprops(cc, 'Centroid', 'BoundingBox');

    % --- Bước 3: Tìm nhóm gần tâm nhất làm gốc ---
    centroids = cat(1, stats.Centroid);
    image_center = [W/2, H/2];
    dist = vecnorm(centroids - image_center, 2, 2);
    [~, idx_center] = min(dist);

    % --- Bước 4: Khởi tạo và gán nhãn ---
    labels = nan(cc.NumObjects, 1);
    labels(idx_center) = 0; % Nhóm gốc đặt nhãn 0

    queue = idx_center; % Hàng đợi để duyệt lan truyền nhãn
    processed_groups = false(cc.NumObjects, 1);
    processed_groups(idx_center) = true;

    % --- Bước 5: Lan truyền nhãn ---
    while ~isempty(queue)
        current_group = queue(1);
        queue(1) = [];

        current_label = labels(current_group);
        pixels = cc.PixelIdxList{current_group};
        [rows, cols] = ind2sub([H, W], pixels);

        visited_gid = []; % Tránh xét lại nhóm cùng vòng lặp

        for i = 1:length(rows)
            r = rows(i);
            c = cols(i);

            % Lan truyền lên trên theo cột
            for y = r-1:-1:1
                gid = labeled_matrix(y, c);
                if gid > 0 && ~processed_groups(gid) && ~ismember(gid, visited_gid)
                    labels(gid) = current_label + 1; % Nhãn tăng dần lên trên
                    queue(end+1) = gid;
                    processed_groups(gid) = true;
                    visited_gid(end+1) = gid;
                    break;
                elseif gid > 0 && processed_groups(gid)
                    break;
                end
            end

            % Lan truyền xuống dưới theo cột
            for y = r+1:H
                gid = labeled_matrix(y, c);
                if gid > 0 && ~processed_groups(gid) && ~ismember(gid, visited_gid)
                    labels(gid) = current_label - 1; % Nhãn giảm dần xuống dưới
                    queue(end+1) = gid;
                    processed_groups(gid) = true;
                    visited_gid(end+1) = gid;
                    break;
                elseif gid > 0 && processed_groups(gid)
                    break;
                end
            end
        end
    end

    % --- Bước 6: Chuẩn hóa nhãn thành số nguyên dương bắt đầu từ 1 ---
    valid_labels = labels(~isnan(labels));

    if isempty(valid_labels)
        warning('Không thể gán nhãn cho bất kỳ vân nào');
        fringe_order = 0;
        fringe_labels = [];
        processed_image = bw_crop;
        return;
    end

    unique_labels = unique(valid_labels);
    labels_new = nan(size(labels));
    for i = 1:length(unique_labels)
        labels_new(labels == unique_labels(i)) = i;
    end
    labels = labels_new;

    % --- Bước 7: Hiển thị kết quả (nếu được yêu cầu) ---
    if display_result
        figure('Name', 'Kết quả gán bậc vân', 'NumberTitle', 'off');
        imshow(bw_crop);
        hold on;

        for k = 1:cc.NumObjects
            if ~isnan(labels(k))
                pixels = cc.PixelIdxList{k};
                [rows, cols] = ind2sub([H, W], pixels);
                coords = [cols, rows]; % [x, y]

                % Tính khoảng cách từ tâm ảnh để đặt nhãn ở vị trí gần tâm nhất
                dists = sqrt((coords(:,1) - image_center(1)).^2 + (coords(:,2) - image_center(2)).^2);
                [~, min_idx] = min(dists);
                label_pos = coords(min_idx, :);

                text(label_pos(1), label_pos(2), num2str(labels(k)), ...
                    'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center');
            end
        end

        title('Gán bậc vân', 'FontSize', 12);
        hold off;
    end

    % --- Bước 8: Trả về kết quả ---
    fringe_order = cc.NumObjects;
    fringe_labels = labels;
    processed_image = bw_crop;

    %     % Hiển thị thống kê
    %     fprintf('Đã phát hiện %d vân\n', fringe_order);
    %     fprintf('Số vân được gán nhãn: %d\n', sum(~isnan(labels)));
    %     if ~isempty(valid_labels)
    %         fprintf('Phạm vi bậc vân: %d đến %d\n', min(unique_labels), max(unique_labels));
    %     end

catch ME
    % Xử lý lỗi
    error_msg = sprintf('Lỗi trong quá trình gán bậc vân:\n%s', ME.message);
    error(error_msg);
end

end
function [recons_surface, figure_handle] = reconSurface_linearPushed(BW, fringe_labels, lambda, tilt_option, show_figure)
% RECONSURFACE_LINEARPUSHED Tái tạo bề mặt 3D từ ảnh vân giao thoa
%
% Cú pháp:
%   [recons_surface, figure_handle] = reconSurface_linearPushed(BW, fringe_labels, lambda, tilt_option, show_figure)
%
% Tham số đầu vào:
%   BW            - Ảnh nhị phân đã cắt biên (logical matrix)
%   fringe_labels - Vector chứa nhãn của các vân (double array)
%   lambda        - Bước sóng ánh sáng (double)
%   tilt_option   - Tùy chọn xử lý ('None', 'Remove tilt', 'Invert', 'Remove + Invert')
%   show_figure   - Có hiển thị figure hay không (logical, optional, default: true)
%
% Tham số đầu ra:
%   recons_surface - Ma trận bề mặt 3D đã tái tạo
%   figure_handle  - Handle của figure (nếu show_figure = true)
%
% Ví dụ:
%   [surface, fig] = reconSurface_linearPushed(BW, [1,2,3,4,5], 632.8e-9, 'Remove tilt');

% Xử lý tham số đầu vào
if nargin < 5
    show_figure = true;
end

% Kiểm tra tham số đầu vào
if isempty(fringe_labels)
    error('Bạn cần gán nhãn vân trước khi nội suy.');
end

if ~islogical(BW)
    error('BW phải là ảnh nhị phân (logical matrix).');
end

% Thiết lập khoảng cách giữa các vân
khoang_cach_van = lambda/2;

% Tìm các thành phần liên thông
cc = bwconncomp(BW);
L = labelmatrix(cc);

% Khởi tạo các mảng điểm 3D
num_labels = max(L(:));
X = []; Y = []; Z = [];

for i = 1:num_labels
    % Tìm các điểm thuộc vân có nhãn i
    [y, x] = find(L == i);

    if i <= length(fringe_labels)
        % Tính độ cao z dựa trên nhãn vân
        z = ones(size(x)) * (fringe_labels(i)) * khoang_cach_van;
        X = [X; x];
        Y = [Y; y];
        Z = [Z; z];
    end
end

% Kiểm tra xem có dữ liệu để nội suy không
if isempty(X)
    error('Không có dữ liệu để nội suy. Kiểm tra lại fringe_labels và BW.');
end

% Nội suy để tạo mặt 3D mượt
[xq, yq] = meshgrid(1:size(BW,2), 1:size(BW,1));
F = scatteredInterpolant(X, Y, Z, 'natural', 'nearest');
Zq = F(xq, yq);
Zq(~isfinite(Zq)) = 0;

% %
% Z_grid_cubic = griddata(X, Y, Z, xq, yq, 'cubic');
% Z_grid_cubic(~isfinite(Z_grid_cubic)) = 0;
%
% % 6. Làm mượt hậu xử lý cho cubic
% Z_cubic_smooth = imgaussfilt(Z_grid_cubic, 2);
% Zq = Z_cubic_smooth;
% %

% Chuyển từ mét sang radian
phi_rad = (4 * pi / lambda) * Zq;
Zq = phi_rad;

% Cắt biên để hiển thị tốt hơn

Z_crop = Zq;


[M, N] = size(Z_crop);
[xGrid, yGrid] = meshgrid(1:N, 1:M);
x = xGrid(:);
y = yGrid(:);
z = Z_crop(:);

% Xử lý theo lựa chọn của người dùng
switch tilt_option
    case 'None'
        Z_processed = Z_crop;

    case 'Remove tilt'
        good = ~isnan(z);
        if sum(good) < 3
            warning('Không đủ điểm hợp lệ để loại bỏ độ nghiêng.');
            Z_processed = Z_crop;
        else
            A = [x, y, ones(size(x))];
            coeff = A(good,:) \ z(good);
            Z_fit = reshape(A * coeff, size(Z_crop));
            Z_processed = Z_crop - Z_fit;
        end

    case 'Invert'
        Z_processed = max(Z_crop(:)) - Z_crop;

    case 'Remove + Invert'
        good = ~isnan(z);
        if sum(good) < 3
            warning('Không đủ điểm hợp lệ để loại bỏ độ nghiêng.');
            Z_leveled = Z_crop;
        else
            A = [x, y, ones(size(x))];
            coeff = A(good,:) \ z(good);
            Z_fit = reshape(A * coeff, size(Z_crop));
            Z_leveled = Z_crop - Z_fit;
        end
        Z_processed = max(Z_leveled(:)) - Z_leveled;

    otherwise
        warning('Tùy chọn không hợp lệ. Sử dụng "None".');
        Z_processed = Z_crop;
end

% Chuẩn hóa bắt đầu từ 0
Z_offset = Z_processed - min(Z_processed(:));

% Gán kết quả đầu ra
recons_surface = Z_offset;

% Hiển thị bề mặt 3D nếu được yêu cầu
if show_figure
    figure_handle = figure;
    surf(xGrid, yGrid, Z_offset);
    shading interp;
    xlabel('X (px)');
    ylabel('Y (px)');
    zlabel('rad');
    title(['3D Surface Linear (Option: ', tilt_option, ')']);
    colormap parula;
    colorbar;
else
    figure_handle = [];
end

end
function [unwrappedPhase, kMap] = unwrapUsingEstimate(estimatedPhase, wrappedPhase)
% Giải Wrapped pha `wrappedPhase` dựa trên pha ước lượng `estimatedPhase`.
% wrappedEstimate = wrapToPi(estimatedPhase);
kMap = round((estimatedPhase - wrappedPhase) / (2*pi));
unwrappedPhase = wrappedPhase + 2*pi * kMap;
% ta có: estimated ~ unwraping_phase
% mà un_phase = wwrapped + k.2pi
% thay 2 vào 1, có: estiamted - wrapped ~ k.2pi
% Suy ra: k ~ (estimated - wrapped)/2pi
end
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
