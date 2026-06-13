function [restored, secret_bits] = ac_extract(marked, side_info, key3)
% AC_EXTRACT  Data extraction and encrypted-share restoration
%             (Sec. IV-C, Fig. 5). Inverse of ac_embed.
%
%   [restored, secret_bits] = ac_extract(marked, side_info, key3)
%
%   marked      : M x N uint8 marked share
%   side_info   : struct from ac_embed (.codes, .lambda_first, .n_secret)
%   key3        : data-hiding key (same seed used in ac_embed)
%
%   restored    : M x N uint8 encrypted share, restored losslessly
%   secret_bits : 1 x n_secret recovered secret bits (PRN-decrypted)
%
%   Per block: type read from the reference-pixel MSB.
%   B1: lambda decoded from the label stream extracted from PREVIOUS
%       blocks (first B1 block's lambda comes from side_info); payload
%       bits harvested; reference pixel restored from the stored MSB;
%       pixels 2-4 restored as ref XOR d (Fig. 5).
%   B2: reference MSB restored from the MSB-information section of the
%       stream in a second pass; pixels 2-4 were never modified.

[M, N] = size(marked);
s = (M * N) / 4;
nbc = N / 2;

% Label code lookup: code string -> effective lambda
code2lam = containers.Map('KeyType', 'char', 'ValueType', 'double');
for l = 0:6
    if ~isempty(side_info.codes{l+1})
        code2lam(side_info.codes{l+1}) = l;
    end
end

buf = zeros(1, 23 * s);      % extracted payload bits (upper bound)
nbuf = 0;
label_ptr = 1;               % next unread bit for label decoding
b1_seen = 0;
b2_list = [];                % B2 block indices in scan order

restored = marked;

for z = 1:s
    br = 2 * floor((z-1) / nbc) + 1;
    bc = 2 * mod(z-1, nbc) + 1;
    ref = double(marked(br, bc));

    if bitget(ref, 8) == 1                      % ── B2 block ──
        b2_list(end+1) = z; %#ok<AGROW>
        continue;
    end

    % ── B1 block ──
    b1_seen = b1_seen + 1;
    if b1_seen == 1
        lam = side_info.lambda_first;
    else
        code = '';
        while true
            if label_ptr > nbuf
                error('ac_extract: label bits not yet available — stream corrupt.');
            end
            code = [code, char('0' + buf(label_ptr))]; %#ok<AGROW>
            label_ptr = label_ptr + 1;
            if isKey(code2lam, code)
                break;
            end
        end
        lam = code2lam(code);
    end
    e1 = lam;
    e2 = 8 - e1;

    pr = [br, br+1, br+1];
    pc = [bc+1, bc, bc+1];
    D = zeros(1, 3);
    ref_msb = 0;
    for u = 1:3
        v = double(marked(pr(u), pc(u)));
        vbits = bitget(v, 8:-1:1);              % MSB first
        if u == 1
            ref_msb = vbits(1);
            payload = vbits(2:e2);              % eps2-1 bits
        else
            payload = vbits(1:e2);              % eps2 bits
        end
        buf(nbuf+1:nbuf+numel(payload)) = payload;
        nbuf = nbuf + numel(payload);
        if e1 > 0
            D(u) = sum(vbits(e2+1:8) .* 2.^(e1-1:-1:0));
        end
    end

    ref_orig = ref_msb * 128 + bitand(ref, 127);
    restored(br, bc) = uint8(ref_orig);
    for u = 1:3
        restored(pr(u), pc(u)) = uint8(bitxor(ref_orig, D(u)));
    end
end

% ── Second pass: restore B2 reference MSBs, then recover secret ─────────────
rest = buf(label_ptr:nbuf);   % stream after the label section

for j = 1:numel(b2_list)
    z = b2_list(j);
    br = 2 * floor((z-1) / nbc) + 1;
    bc = 2 * mod(z-1, nbc) + 1;
    msb = rest(j);
    restored(br, bc) = uint8(msb * 128 + bitand(double(marked(br, bc)), 127));
end

n_sec = side_info.n_secret;
sec_enc = rest(numel(b2_list)+1 : numel(b2_list)+n_sec);

% Decrypt with PRN(key3)
rng(key3, 'twister');
prn = randi([0, 1], 1, n_sec);
secret_bits = bitxor(sec_enc, prn);
end
