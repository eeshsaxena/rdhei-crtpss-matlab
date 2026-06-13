function img = crtp_reconstruct(shares_subset, q_subset, q0, key1, s_total)
% CRTP_RECONSTRUCT  Recover the (scrambled) image from any k encrypted
%                   shares via CRT for polynomials (Eq. 14).
%
%   img = crtp_reconstruct(shares_subset, q_subset, q0, key1, s_total)
%
%   shares_subset : 1 x k cell of M x N uint8 restored encrypted shares
%   q_subset      : 1 x k share moduli matching the shares (uint64)
%   q0            : auxiliary modulus (uint64)
%   key1          : seed used for R1 during sharing
%   s_total       : total number of 2x2 blocks (regenerates R1 sequence)
%
%   Per pixel:
%     Y(x)    = XOR_i [ m^e_i(x) (*) Q_i(x) (*) Q_i^{-1}(x) ]  mod Q(x)
%     m_zu(x) = [ Y(x) mod q0(x) ]  XOR  R1_z
%
%   where Q = prod(q_i), Q_i = Q/q_i, Q_i^{-1} = inverse of Q_i mod q_i.
%   R2 vanishes automatically: (R2 (*) q0) mod q0 = 0.

k = numel(shares_subset);
[M, N] = size(shares_subset{1});

% Precompute CRT constants
Q = uint64(1);
for i = 1:k
    Q = gf2_mul(Q, q_subset(i));
end
Qi  = cell(1, k);
inv = cell(1, k);
for i = 1:k
    [Qi{i}, rem] = gf2_divmod(Q, q_subset(i));
    if rem ~= 0
        error('q_subset(%d) does not divide Q — moduli inconsistent.', i);
    end
    inv{i} = gf2_inv(gf2_mod(Qi{i}, q_subset(i)), q_subset(i));
    % Pre-multiply Qi * inv (degree < deg(Q)+8, reduced per-pixel below)
    Qi{i} = gf2_mul(Qi{i}, inv{i});
end

% Regenerate R1 from key1
rng(key1, 'twister');
R1 = uint64(randi([0, 255], 1, s_total));

nbc = N / 2;
offr = [0 0 1 1];
offc = [0 1 0 1];

img = zeros(M, N, 'uint8');
for z = 1:s_total
    br = 2 * floor((z-1) / nbc) + 1;
    bc = 2 * mod(z-1, nbc) + 1;
    for u = 1:4
        Y = uint64(0);
        for i = 1:k
            me = uint64(shares_subset{i}(br + offr(u), bc + offc(u)));
            Y = bitxor(Y, gf2_mod(gf2_mul(me, Qi{i}), Q));
        end
        m = bitxor(gf2_mod(Y, q0), R1(z));
        img(br + offr(u), bc + offc(u)) = uint8(m);
    end
end
end
