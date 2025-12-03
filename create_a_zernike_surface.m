% my_zernike
% hard-code 36 he so dau Zernike

function [surface, pos_coeff, val_coeff] = create_a_zernike_surface()
gridsize = 512;
gridsize = gridsize*sqrt(2);
x = linspace(-1,1, gridsize);
y = linspace(-1,1, gridsize);
[X, Y] = meshgrid(x,y);

[theta, rho] = cart2pol(X,Y);
% TẠO MASK: Chỉ lấy giá trị trong hình tròn bán kính 1
mask = rho <= 1;

% Mảng 1: Vị trí (Index) các Zernike 
% pos_coeff = [1, 2, 3, 4, 5, 6, 7 ]; 
% val_coeff = [];

pos_coeff = 4:16;      % index các đa thức Zernike
% max_amp = 3;
% val_coeff = (rand(1, length(pos_coeff)) - 0.5) * 2 * max_amp;
val_coeff =[ ...
     3.0, ...   % 4   (PEAK lồi mạnh)
    -2.50, ...   % 5   (LÕM sâu)
     0.90, ...   % 6
     2.20, ...   % 7   (lồi mạnh)
    -1.70, ...   % 8   (lõm)
     1.40, ...   % 9
    -2.20, ...   % 10  (lõm mạnh)
     0.85, ...   % 11
     1.10, ...   % 12
    -1.90, ...   % 13  (lõm rõ)
     0.70, ...   % 14
    -1.30, ...   % 15
     2.00  ...   % 16  (lồi mạnh)
];

% Kiểm tra an toàn: Hai mảng phải có cùng độ dài
if length(pos_coeff) ~= length(val_coeff)
    error('Lỗi: Số lượng vị trí và số lượng hệ số không bằng nhau!');
end
surface = zeros(size(rho));
for k = 1:length(pos_coeff)
    j = pos_coeff(k); 
    c = val_coeff(k);  
    if c==0
        continue;
    end
    Z_term = my_get_zernike_poly(j, rho, theta);
    surface = surface + c*Z_term;
end

x_square = abs(X) <= sqrt(2)/2;
y_square = abs(Y) <= sqrt(2)/2;
mask_square = x_square & y_square;
surface_square = surface;   % hình vuông từ dữ liệu gốc (không NaN)
surface_square(~mask_square) = NaN;
% Tìm tất cả vị trí không phải NaN
[row_idx, col_idx] = find(~isnan(surface_square));

% Lấy bounding box của vùng có dữ liệu
rmin = min(row_idx);
rmax = max(row_idx);
cmin = min(col_idx);
cmax = max(col_idx);

% Cắt ra vùng vuông sạch, không NaN
surface_square = surface_square(rmin:rmax, cmin:cmax);
surface = surface_square;
end