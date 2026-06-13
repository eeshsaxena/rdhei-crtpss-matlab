%% MAIN_DEMO  Full pipeline demo: RDH with secret encrypted image sharing
%             and adaptive coding (Fang et al., IEEE IoTJ 2025).
%
%  Pipeline (Fig. 1):
%   Content owner : scramble (Key0) -> CRTP-SS (k,n)-threshold -> n shares
%   Data hider i  : adaptive-coding embed (Key3_i) -> marked share i
%   Receiver      : any k marked shares -> extract + restore shares ->
%                   CRT-P reconstruct -> unscramble -> original image
%
%  Uses the paper's (3,4)-threshold and moduli
%  q0=111110101, q1=100011011, q2=100011101, q3=100101011, q4=100101101.

clear; clc; close all;
addpath(genpath('src'));
addpath(genpath('utils'));

%% ── Parameters ──────────────────────────────────────────────────────────────
k = 3;  n = 4;                       % (3,4)-threshold, as in Sec. V-A
M = 64; N = 64;                      % image size (use 512x512 for paper scale)
key0 = [11, 22];                     % [block scramble, pixel scramble]
key1 = 33;                           % R1 sequence
key2 = 44;                           % R2 sequence
key3 = [101, 102, 103, 104];         % per-data-hider embedding keys

[q0, q_list] = irreducible_polys(n);

rng(2025);
img = make_test_gray_image('smooth', M, N);

%% ── Content owner: scramble + CRTP-SS sharing ───────────────────────────────
fprintf('=== Content owner: scrambling + CRTP-SS (%d,%d)-threshold ===\n', k, n);
t0 = tic;
scr = scramble_image(img, key0);
shares = crtp_share(scr, n, k, key1, key2, q_list, q0);
fprintf('    Encryption time: %.1f ms\n', toc(t0)*1000);
fprintf('    Original entropy: %.4f | Share entropies:', ...
        image_metrics('entropy', img));
for i = 1:n
    fprintf(' %.4f', image_metrics('entropy', shares{i}));
end
fprintf('  (ideal: 8)\n');

%% ── Data hiders: adaptive-coding embedding ──────────────────────────────────
fprintf('\n=== Data hiders: adaptive-coding embedding ===\n');

% First query capacity with an empty payload, then fill ~90%% of pure EC
[~, ~, stats0] = ac_embed(shares{1}, [], key3(1));
fprintf('    ECs=%d, label=%d, MSB=%d, pure EC=%d bits, pure ER=%.4f bpp\n', ...
        stats0.ECs, stats0.label_bits, stats0.msb_bits, ...
        stats0.pureEC, stats0.pureER);

n_secret = floor(stats0.pureEC * 0.9);
secrets = cell(1, n);
marked  = cell(1, n);
sides   = cell(1, n);
for i = 1:n
    secrets{i} = randi([0, 1], 1, n_secret);
    [marked{i}, sides{i}, st] = ac_embed(shares{i}, secrets{i}, key3(i));
    fprintf('    Hider %d: embedded %d secret bits (pure ER %.4f bpp)\n', ...
            i, n_secret, st.pureER);
end

%% ── Receiver: any k=3 of n=4 marked shares ──────────────────────────────────
fprintf('\n=== Receiver: extraction + recovery from k=%d shares ===\n', k);
s_total = (M * N) / 4;
combos = nchoosek(1:n, k);
all_ok = true;
for c = 1:size(combos, 1)
    idx = combos(c, :);
    restored = cell(1, k);
    sec_ok = true;
    for j = 1:k
        [restored{j}, sec] = ac_extract(marked{idx(j)}, sides{idx(j)}, key3(idx(j)));
        sec_ok = sec_ok && isequal(sec, secrets{idx(j)});
    end
    rec_scr = crtp_reconstruct(restored, q_list(idx), q0, key1, s_total);
    rec = unscramble_image(rec_scr, key0);
    lossless = isequal(rec, img);
    all_ok = all_ok && lossless && sec_ok;
    fprintf('    Shares %s | secrets exact: %d | image lossless: %d | PSNR: %s\n', ...
            mat2str(idx), sec_ok, lossless, ...
            num2str(image_metrics('psnr', img, rec)));
end
if all_ok
    fprintf('    >>> ALL %d combinations: perfect recovery <<<\n', size(combos,1));
end

%% ── Visualisation ───────────────────────────────────────────────────────────
idx = [1 2 3];
restored = cell(1, k);
for j = 1:k
    [restored{j}, ~] = ac_extract(marked{idx(j)}, sides{idx(j)}, key3(idx(j)));
end
rec = unscramble_image(crtp_reconstruct(restored, q_list(idx), q0, key1, s_total), key0);

figure('Name', 'RDHEI-CRTPSS Pipeline', 'NumberTitle', 'off', ...
       'Position', [60 60 1280 640]);
subplot(2, 4, 1); imshow(img);        title('Original');
subplot(2, 4, 2); imshow(scr);        title('Scrambled (Key_0)');
subplot(2, 4, 3); imshow(shares{1});  title('Encrypted share 1');
subplot(2, 4, 4); imshow(marked{1});  title('Marked share 1');
subplot(2, 4, 5); imshow(marked{2});  title('Marked share 2');
subplot(2, 4, 6); imshow(marked{3});  title('Marked share 3');
subplot(2, 4, 7); imshow(rec);        title('Recovered (lossless)');
subplot(2, 4, 8);
histogram(double(shares{1}(:)), 0:255); title('Share 1 histogram (uniform)');
sgtitle(sprintf('RDH with CRTP-SS + Adaptive Coding — (%d,%d)-threshold', k, n));
