%% RUN_EXPERIMENTS  Reproduce the paper's experimental analyses (Sec. V):
%   E1: Embedding statistics per image type (Table II equivalent:
%       ECs, label info, MSB info, pure EC, pure ER)
%   E2: Security — information entropy of shares (Table VI equivalent)
%   E3: Security — key sensitivity via NPCR / UACI (Table VII equivalent)
%   E4: Histograms of original vs. encrypted shares (Fig. 10 equivalent)
%
%  Synthetic images are used; swap in imread('house.png') etc. to run on
%  the paper's six 512x512 test images.

clear; clc; close all;
addpath(genpath('src'));
addpath(genpath('utils'));

k = 3;  n = 4;
M = 128; N = 128;
[q0, q_list] = irreducible_polys(n);
key0 = [11, 22]; key1 = 33; key2 = 44; key3 = 101;

types = {'smooth', 'mixed', 'complex'};

%% ── E1: Embedding statistics ────────────────────────────────────────────────
fprintf('=== E1: Embedding statistics (Table II equivalent) ===\n');
fprintf('%-10s %10s %12s %10s %10s %10s\n', ...
        'Image', 'ECs(bits)', 'Label(bits)', 'MSB(bits)', 'PureEC', 'PureER');
fprintf('%s\n', repmat('-', 1, 68));

share1_by_type = cell(1, numel(types));
img_by_type    = cell(1, numel(types));
for t = 1:numel(types)
    rng(t * 11);
    img = make_test_gray_image(types{t}, M, N);
    scr = scramble_image(img, key0);
    shares = crtp_share(scr, n, k, key1, key2, q_list, q0);
    [~, ~, st] = ac_embed(shares{1}, [], key3);
    fprintf('%-10s %10d %12d %10d %10d %9.4f\n', types{t}, ...
            st.ECs, st.label_bits, st.msb_bits, st.pureEC, st.pureER);
    share1_by_type{t} = shares;
    img_by_type{t}    = img;
end

%% ── E2: Information entropy ─────────────────────────────────────────────────
fprintf('\n=== E2: Information entropy (Table VI equivalent, ideal = 8) ===\n');
fprintf('%-10s %10s %8s %8s %8s %8s %9s\n', ...
        'Image', 'Original', 'I1e', 'I2e', 'I3e', 'I4e', 'Average');
fprintf('%s\n', repmat('-', 1, 66));
for t = 1:numel(types)
    img = img_by_type{t};
    shares = share1_by_type{t};
    e = zeros(1, n);
    for i = 1:n
        e(i) = image_metrics('entropy', shares{i});
    end
    fprintf('%-10s %10.4f %8.4f %8.4f %8.4f %8.4f %9.4f\n', ...
            types{t}, image_metrics('entropy', img), e, mean(e));
end

%% ── E3: Key sensitivity — NPCR / UACI ───────────────────────────────────────
fprintf('\n=== E3: Key sensitivity, Key1 -> Key1+1 (Table VII equivalent) ===\n');
fprintf('    Ideal: NPCR = 99.6094%%, UACI = 33.4635%%\n');
fprintf('%-10s %8s %12s %12s\n', 'Image', 'Share', 'NPCR(%)', 'UACI(%)');
fprintf('%s\n', repmat('-', 1, 46));
for t = 1:numel(types)
    img = img_by_type{t};
    scr = scramble_image(img, key0);
    shA = crtp_share(scr, n, k, key1,     key2, q_list, q0);
    shB = crtp_share(scr, n, k, key1 + 1, key2, q_list, q0);   % 1-key change
    npcr_all = zeros(1, n); uaci_all = zeros(1, n);
    for i = 1:n
        npcr_all(i) = image_metrics('npcr', shA{i}, shB{i});
        uaci_all(i) = image_metrics('uaci', shA{i}, shB{i});
        fprintf('%-10s %8d %12.4f %12.4f\n', types{t}, i, npcr_all(i), uaci_all(i));
    end
    fprintf('%-10s %8s %12.4f %12.4f   (averages)\n', '', 'avg', ...
            mean(npcr_all), mean(uaci_all));
end

%% ── E4: Histograms ──────────────────────────────────────────────────────────
figure('Name', 'Histograms (Fig. 10 equivalent)', 'NumberTitle', 'off', ...
       'Position', [80 80 1280 320]);
img = img_by_type{1};
shares = share1_by_type{1};
subplot(1, n+1, 1);
histogram(double(img(:)), 0:255); title('Original (smooth)');
xlabel('Value'); ylabel('Frequency');
for i = 1:n
    subplot(1, n+1, 1+i);
    histogram(double(shares{i}(:)), 0:255);
    title(sprintf('Share %d', i));
    xlabel('Value');
end
sgtitle('Pixel histograms: original vs. CRTP-SS encrypted shares (uniform 0-255)');

fprintf('\nDone. Encrypted-share histograms should be uniform over 0-255,\n');
fprintf('unlike CRTI-SS where pixels cannot reach values above the moduli.\n');
