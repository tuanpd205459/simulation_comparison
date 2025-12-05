clc; clear; close all;

%% ===================== 1. Khởi tạo tham số =======================

N = 300;                    % kích thước lưới
max_amp = 3;                % biên độ lồi - lõm (càng lớn càng gồ ghề)
idx = 4:17;                 % chọn mode Zernike để tạo lỗi

% tạo hệ số bao gồm cả dương + âm => chắc chắn có lồi và lõm
n = length(idx);
half = floor(n/2);
rng('shuffle');             % random mỗi lần chạy khác nhau

pos =  (rand(1,half) * max_amp);
neg = -(rand(1,n-half) * max_amp);
val_coeff = [pos neg];
val_coeff = val_coeff(randperm(n));  % xáo trộn ngẫu nhiên

disp('Hệ số Zernike sử dụng:');
disp(val_coeff);

%% ===================== 2. Tạo lưới cực và mặt nạ =======================

[x,y] = meshgrid(linspace(-1,1,N));
[r,theta] = cart2pol(x,y);

mask = r <= 1;              % mặt nạ pupil tròn
r_norm = r;                 % r đã chuẩn hóa sẵn vì [-1..1]

%% ===================== 3. Tính bề mặt Zernike =======================

Z = zeros(N);
k = 1;

for n_id = idx
    Zk = zernikeMode(n_id, r_norm, theta);
    Z = Z + val_coeff(k) * Zk .* mask;
    k = k + 1;
end

%% ===================== 4. Hiển thị kết quả =======================

figure; 
subplot(1,2,1);
imagesc(Z); axis image; colormap jet; colorbar;
title('Surface Zernike (2D)');

subplot(1,2,2);
surf(Z, 'EdgeColor','none'); colormap jet; colorbar;
title('Surface Zernike (3D)');
view(45,35)
camlight; lighting phong;

function Z = zernikeMode(j, r, theta)
% j: Zernike index (Noll indexing)
% r,theta dạng meshgrid

[n,m] = noll_to_nm(j);    % chuyển từ Noll → (n,m)

if m==0
    R = radialZernike(n,0,r);
    Z = R;
elseif m>0
    R = radialZernike(n,m,r);
    Z = R.*cos(m*theta);
else
    R = radialZernike(n,-m,r);
    Z = R.*sin(-m*theta);
end
Z(r>1)=0;
end
function R = radialZernike(n,m,r)
R = zeros(size(r));
for k=0:(n-m)/2
    R = R + (-1)^k * factorial(n-k)/( factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k) ) ...
        .* r.^(n-2*k);
end
end
function [n,m] = noll_to_nm(j)
n = floor((sqrt(8*j-7)-1)/2);
p = j - n.*(n+1)/2;
m = 2*p - 1 - mod(n,2);
end
