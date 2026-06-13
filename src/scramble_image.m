function scr = scramble_image(img, key0)
% SCRAMBLE_IMAGE  Block + in-block pixel scrambling with Key0 (Sec. IV-A).
%
%   scr = scramble_image(img, key0)
%
%   img  : M x N uint8 grayscale image, M and N even (2x2 blocks)
%   key0 : 1x2 vector [Key0_1, Key0_2]
%          Key0_1 seeds the block permutation, Key0_2 the per-block
%          pixel permutations.
%
%   Block z (row-major over 2x2 tiles) is moved to position perm(z) and
%   its 4 pixels [tl tr bl br] are permuted by P(z,:).
%   Inverse operation: unscramble_image.

[M, N] = size(img);
if mod(M,2) || mod(N,2)
    error('Image dimensions must be even for 2x2 blocks.');
end
s = (M * N) / 4;

rng(key0(1), 'twister');
perm = randperm(s);                 % source block z -> destination perm(z)
rng(key0(2), 'twister');
P = zeros(s, 4);
for z = 1:s
    P(z, :) = randperm(4);          % pixel permutation of source block z
end

nbc = N / 2;                        % blocks per row
% pixel offsets within a block: 1=tl, 2=tr, 3=bl, 4=br
offr = [0 0 1 1];
offc = [0 1 0 1];

scr = zeros(M, N, 'uint8');
for z = 1:s
    sr = 2 * floor((z-1) / nbc) + 1;       % source block top-left
    sc = 2 * mod(z-1, nbc) + 1;
    dz = perm(z);
    dr = 2 * floor((dz-1) / nbc) + 1;      % destination block top-left
    dc = 2 * mod(dz-1, nbc) + 1;
    for u = 1:4
        scr(dr + offr(u), dc + offc(u)) = ...
            img(sr + offr(P(z,u)), sc + offc(P(z,u)));
    end
end
end
