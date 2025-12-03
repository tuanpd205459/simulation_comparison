function ph_L2 = unwrap_phase_LS(unph, value_x, value_y)
% UNWRAP_PHASE_L2 - Giải pha từ gradient đã hiệu chỉnh (2D)
%
% Syntax:
%   ph_L2 = unwrap_phase_L2(unph, value_x, value_y)
%
% Inputs:
%   unph     - Pha gói (2D array)
%   value_x  - Kích thước hoặc thông số x cho psf2otf_test
%   value_y  - Kích thước hoặc thông số y cho psf2otf_test
%
% Output:
%   ph_L2    - Pha đã được giải (unwrapped)

    % --- 1. Tạo toán tử vi phân trong miền Fourier ---
    dx = psf2otf_test([-1,1;0,0],[value_x,value_y]); % Đạo hàm theo x
    dy = psf2otf_test([-1,0;1,0],[value_x,value_y]); % Đạo hàm theo y

    % --- 2. Tính tổng bình phương đạo hàm ---
    DTD = abs(dx).^2 + abs(dy).^2;

    % --- 3. Gradient pha gói ---
    dadx = real(ifft2(fft2(unph).*dx));
    dady = real(ifft2(fft2(unph).*dy));

    % --- 4. Loại bỏ bước nhảy 2*pi trong gradient ---
    dadx_G = dadx - pi*round(dadx/pi);
    dady_G = dady - pi*round(dady/pi);

    % --- 5. Phục hồi pha liên tục ---
    ph_L2 = real(ifft2((fft2(dadx_G).*conj(dx) + fft2(dady_G).*conj(dy)) ./ (DTD + eps)));

end
