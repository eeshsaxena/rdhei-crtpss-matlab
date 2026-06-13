function img = make_test_gray_image(type, M, N)
% MAKE_TEST_GRAY_IMAGE  Synthetic 8-bit grayscale test images.
%
%   img = make_test_gray_image(type, M, N)
%
%   type : 'smooth'   - gradient + low-frequency sinusoid (mostly B1 blocks,
%                       high embedding rate, like House/Airplane)
%          'complex'  - uniform random noise (mostly B2 blocks, low rate,
%                       like Baboon texture)
%          'mixed'    - smooth left half, noisy right half
%
%   Replace with imread(...) for real test images (House, Baboon, ...).

[cc, rr] = meshgrid(1:N, 1:M);
switch lower(type)
    case 'smooth'
        img = uint8(mod(round(80 + 60*sin(rr/17) + 50*cos(cc/23) + ...
                              (rr + cc)/8), 256));
    case 'complex'
        img = uint8(randi([0, 255], M, N));
    case 'mixed'
        img = uint8(mod(round(80 + 60*sin(rr/17) + (rr + cc)/8), 256));
        noise = uint8(randi([0, 255], M, N));
        img(:, floor(N/2)+1:end) = noise(:, floor(N/2)+1:end);
    otherwise
        error('Unknown type: %s', type);
end
end
