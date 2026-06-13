function img = unscramble_image(scr, key0)
% UNSCRAMBLE_IMAGE  Inverse of scramble_image (Key0 re-scrambling, Sec. IV-C).
%
%   img = unscramble_image(scr, key0)
%
%   Regenerates the same block permutation and per-block pixel permutations
%   from key0 and applies them in reverse.

[M, N] = size(scr);
s = (M * N) / 4;

rng(key0(1), 'twister');
perm = randperm(s);
rng(key0(2), 'twister');
P = zeros(s, 4);
for z = 1:s
    P(z, :) = randperm(4);
end

nbc = N / 2;
offr = [0 0 1 1];
offc = [0 1 0 1];

img = zeros(M, N, 'uint8');
for z = 1:s
    sr = 2 * floor((z-1) / nbc) + 1;
    sc = 2 * mod(z-1, nbc) + 1;
    dz = perm(z);
    dr = 2 * floor((dz-1) / nbc) + 1;
    dc = 2 * mod(dz-1, nbc) + 1;
    for u = 1:4
        img(sr + offr(P(z,u)), sc + offc(P(z,u))) = ...
            scr(dr + offr(u), dc + offc(u));
    end
end
end
