function [marked, side_info, stats] = ac_embed(share, secret_bits, key3)
% AC_EMBED  Adaptive-coding data embedding into an encrypted share
%           (Sec. IV-B: XORP + Huffman label coding, Fig. 4).
%
%   [marked, side_info, stats] = ac_embed(share, secret_bits, key3)
%
%   share       : M x N uint8 encrypted share (2x2 blocks)
%   secret_bits : row vector of 0/1 secret bits (encrypted internally
%                 with a PRN sequence seeded by key3)
%   key3        : data-hiding key (scalar seed)
%
%   marked      : M x N uint8 marked share
%   side_info   : struct passed to ac_extract:
%                   .codes        1x7 cell, Huffman label code per lambda
%                   .lambda_first effective lambda of the first B1 block
%                   .n_secret     number of secret bits
%   stats       : struct with .ECs, .label_bits, .msb_bits, .pureEC,
%                 .pureER (bpp) — mirrors Table II / Eq. 15.
%
%   Block classification (XORP, Eq. 13): d_zu = m^e_z1 XOR m^e_zu,
%   Dmax = max(d). Dmax <= 63 -> B1 (embeddable, lambda by Table I),
%   else B2 (nonembeddable). The reference-pixel MSB marks the type
%   (0 = B1, 1 = B2).
%
%   B1 block bit layout (lambda -> eps1 = lambda, eps2 = 8 - eps1):
%     ref     : [type 0][7 original LSBs]
%     pixel 2 : [original ref MSB][eps2-1 payload][eps1 bits of d_z2]
%     pixel 3 : [eps2 payload][eps1 bits of d_z3]
%     pixel 4 : [eps2 payload][eps1 bits of d_z4]
%   Capacity per B1 block: 3*eps2 - 1 bits.
%
%   Payload stream = [label codes of B1 blocks 2..end]
%                  + [original MSBs of B2 reference pixels]
%                  + [encrypted secret bits], zero-padded to capacity.
%   Decodability: max label length <= 5 <= min block capacity, so the
%   receiver always holds block j's label before reaching block j.

[M, N] = size(share);
s = (M * N) / 4;
nbc = N / 2;
offr = [0 0 1 1];
offc = [0 1 0 1];

% ── Pass 1: classify blocks ─────────────────────────────────────────────────
typ  = zeros(1, s);          % 1 = B1, 2 = B2
lam  = -ones(1, s);          % raw lambda for B1 blocks
Dall = zeros(s, 3);
refs = zeros(1, s);
lam_hist = zeros(1, 7);

for z = 1:s
    br = 2 * floor((z-1) / nbc) + 1;
    bc = 2 * mod(z-1, nbc) + 1;
    p = double([share(br,bc), share(br,bc+1), share(br+1,bc), share(br+1,bc+1)]);
    refs(z) = p(1);
    D = bitxor(p(1), p(2:4));
    Dmax = max(D);
    Dall(z, :) = D;
    if Dmax <= 63
        typ(z) = 1;
        if Dmax == 0
            lam(z) = 0;
        else
            lam(z) = floor(log2(Dmax)) + 1;   % ceil(log2(Dmax+1))
        end
        lam_hist(lam(z)+1) = lam_hist(lam(z)+1) + 1;
    else
        typ(z) = 2;
    end
end

% ── Huffman label codes with merge rule ─────────────────────────────────────
[codes, merge_map] = huffman_label(lam_hist);

b1_idx = find(typ == 1);
lam_eff = -ones(1, s);
for z = b1_idx
    lam_eff(z) = merge_map(lam(z)+1);
end

% ── Build payload stream ────────────────────────────────────────────────────
label_stream = [];
for j = 2:numel(b1_idx)
    z = b1_idx(j);
    label_stream = [label_stream, codes{lam_eff(z)+1} - '0']; %#ok<AGROW>
end

b2_idx = find(typ == 2);
msb_stream = zeros(1, numel(b2_idx));
for j = 1:numel(b2_idx)
    msb_stream(j) = bitget(refs(b2_idx(j)), 8);
end

% Encrypt secret bits with PRN(key3)
secret_bits = double(secret_bits(:)');
rng(key3, 'twister');
prn = randi([0, 1], 1, numel(secret_bits));
sec_enc = bitxor(secret_bits, prn);

caps = 3 * (8 - lam_eff(b1_idx)) - 1;
ECs  = sum(caps);
stream = [label_stream, msb_stream, sec_enc];

if numel(stream) > ECs
    error(['ac_embed: payload (%d bits) exceeds capacity (%d bits). ', ...
           'Pure EC available for secret data: %d bits.'], ...
          numel(stream), ECs, ECs - numel(label_stream) - numel(msb_stream));
end
stream = [stream, zeros(1, ECs - numel(stream))];   % zero-pad

% ── Pass 2: write marked blocks ─────────────────────────────────────────────
marked = share;
pos = 1;
pw = 2.^(7:-1:0);
for z = 1:s
    br = 2 * floor((z-1) / nbc) + 1;
    bc = 2 * mod(z-1, nbc) + 1;
    if typ(z) == 2
        % B2: set type bit 1, pixels 2-4 untouched
        marked(br, bc) = uint8(bitor(bitand(refs(z), 127), 128));
        continue;
    end
    e1 = lam_eff(z);
    e2 = 8 - e1;
    marked(br, bc) = uint8(bitand(refs(z), 127));   % type bit 0
    pr = [br, br+1, br+1];                          % pixel 2,3,4 rows
    pc = [bc+1, bc, bc+1];                          % pixel 2,3,4 cols
    for u = 1:3
        if u == 1
            head = [bitget(refs(z), 8), stream(pos:pos+e2-2)];
            pos = pos + e2 - 1;
        else
            head = stream(pos:pos+e2-1);
            pos = pos + e2;
        end
        if e1 > 0
            dbits = bitget(Dall(z,u), e1:-1:1);
        else
            dbits = [];
        end
        marked(pr(u), pc(u)) = uint8(sum([head, dbits] .* pw));
    end
end

% ── Side info and stats ─────────────────────────────────────────────────────
if isempty(b1_idx)
    side_info.lambda_first = -1;
else
    side_info.lambda_first = lam_eff(b1_idx(1));
end
side_info.codes    = codes;
side_info.n_secret = numel(secret_bits);

stats.ECs        = ECs;
stats.label_bits = numel(label_stream);
stats.msb_bits   = numel(msb_stream);
stats.pureEC     = ECs - stats.label_bits - stats.msb_bits;
stats.pureER     = stats.pureEC / (M * N);
stats.n_B1       = numel(b1_idx);
stats.n_B2       = numel(b2_idx);
end
