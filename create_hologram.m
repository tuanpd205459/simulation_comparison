function [hologram, wrapped_phase, carrier] = create_hologram(surface, noise_level)
N = size(surface,1);
M = size(surface,2);

object_phase_without_noise = surface;

% Thêm nhiễu vào pha đối tượng
sigma_signal = std(object_phase_without_noise(:), 'omitnan');

% 3. Tạo nhiễu (Bỏ 'omitnan' trong randn)
noise = (noise_level .* sigma_signal) .* randn(size(object_phase_without_noise));

object_phase = object_phase_without_noise + noise;
%% 3. TẠO HOLOGRAM
fx = 30 / N; % Tần số sóng mang
fy = -30 / M;

[X, Y] = meshgrid(1:N, 1:M);

% Cường độ nền và điều biến
a = 1.0; % Background intensity
b = 0.8; % Modulation depth

% Sóng mang phẳng (plane wave carrier)
carrier = 2 * pi * (fx * X + fy * Y);

% --- Tạo hologram (ảnh giao thoa) theo công thức mới ---
hologram = a + b .* cos(carrier + object_phase);
hologram = mat2gray(hologram);

wrapped_phase = wrapToPi(object_phase);

end