function shares = crtp_share(img, n, k, key1, key2, q_list, q0)
% CRTP_SHARE  CRTP-SS encryption: split image into n encrypted shares
%             (Sec. III-B and IV-A, Eqs. 5-6).
%
%   shares = crtp_share(img, n, k, key1, key2, q_list, q0)
%
%   img    : M x N uint8 (already scrambled), M and N even
%   n, k   : (k,n)-threshold; any k shares reconstruct
%   key1   : seed for per-block 8-bit random sequence R1
%   key2   : seed for per-block random sequence R2
%   q_list : 1 x n share moduli (irreducible degree-8 polys, uint64)
%   q0     : auxiliary modulus (irreducible degree-8 poly, uint64)
%
%   shares : 1 x n cell of M x N uint8 encrypted shares
%
%   Per pixel m_zu of block z (Eq. 5):
%     Y(x)        = m_zu(x) XOR R1_z XOR R2_z (*) q0(x)
%     m^e_zui(x)  = Y(x) mod q_i(x)                       (Eq. 6)
%
%   R2 bit-width: deg(Y) must be < 8k for unique CRT reconstruction,
%   so width = min(16, 8*(k-1)). The paper's (3,4) setting gives 16 bits.
%
%   XOR-preservation (Eq. 8): within a block,
%     m^e_z1i XOR m^e_zui == m_z1 XOR m_zu
%   which vacates room for the adaptive-coding embedding.

[M, N] = size(img);
s = (M * N) / 4;

w2 = min(16, 8 * (k - 1));          % R2 bit width (16 for k>=3, 8 for k=2)
rng(key1, 'twister');
R1 = uint64(randi([0, 255], 1, s));
rng(key2, 'twister');
R2 = uint64(randi([0, 2^w2 - 1], 1, s));

nbc = N / 2;
offr = [0 0 1 1];
offc = [0 1 0 1];

shares = cell(1, n);
for i = 1:n
    shares{i} = zeros(M, N, 'uint8');
end

for z = 1:s
    br = 2 * floor((z-1) / nbc) + 1;
    bc = 2 * mod(z-1, nbc) + 1;
    r2q0 = gf2_mul(R2(z), q0);      % R2_z (*) q0(x), same for all 4 pixels
    for u = 1:4
        m = uint64(img(br + offr(u), bc + offc(u)));
        Y = bitxor(bitxor(m, R1(z)), r2q0);
        for i = 1:n
            shares{i}(br + offr(u), bc + offc(u)) = ...
                uint8(gf2_mod(Y, q_list(i)));
        end
    end
end
end
