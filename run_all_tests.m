%% RUN_ALL_TESTS  Unit tests for the CRTP-SS + adaptive-coding RDHEI scheme.
%
%  T1: GF(2)[x] arithmetic (mul/divmod identity, modular inverse)
%  T2: Paper's worked example (Sec. IV-A, Eqs. 9-12):
%      block [192 191; 195 180], R1=128, R2=1, q0=111110101, q1=100011011
%      must give encrypted share-1 block [174 209; 173 218]
%  T3: XOR-preservation property (Eq. 8)
%  T4: CRT-P pixel reconstruction from every k-subset of shares
%  T5: scramble/unscramble inverse
%  T6: AC embed/extract round-trip (share restored + secret bit-exact)
%  T7: full pipeline lossless for all C(4,3) combos, smooth & mixed images
%  T8: complex (B2-heavy) image still restores losslessly
%  T9: capacity-overflow error is raised

clear; clc;
addpath(genpath('src'));
addpath(genpath('utils'));

results = cell(0, 2);

%% ── T1: GF(2)[x] arithmetic ─────────────────────────────────────────────────
fprintf('--- T1: GF(2)[x] arithmetic ---\n');
q_aes = uint64(bin2dec('100011011'));
ok = true;
rng(0);
for a = [1, 2, 3, 83, 202, 255]
    inv = gf2_inv(uint64(a), q_aes);
    ok = ok && (gf2_mod(gf2_mul(uint64(a), inv), q_aes) == 1);
end
results(end+1,:) = {ok, 'T1: modular inverse in GF(2^8)'};

ok = true;
for t = 1:200
    a = uint64(randi([1, 2^20]));
    b = uint64(randi([1, 2^10]));
    [q, r] = gf2_divmod(a, b);
    ok = ok && (bitxor(gf2_mul(q, b), r) == a) && ...
         (r == 0 || gf2_deg(r) < gf2_deg(b));
end
results(end+1,:) = {ok, 'T1: divmod identity a = q*b + r'};

%% ── T2: paper worked example ────────────────────────────────────────────────
fprintf('--- T2: paper worked example ---\n');
[q0, q_list] = irreducible_polys(4);
R1 = uint64(128);  R2 = uint64(1);
block    = [192, 191, 195, 180];
expected = [174, 209, 173, 218];
got = zeros(1, 4);
for u = 1:4
    Y = bitxor(bitxor(uint64(block(u)), R1), gf2_mul(R2, q0));
    got(u) = double(gf2_mod(Y, q_list(1)));
end
results(end+1,:) = {isequal(got, expected), ...
    sprintf('T2: block [192 191 195 180] -> share1 %s', mat2str(got))};

%% ── T3: XOR preservation ────────────────────────────────────────────────────
ok = true;
for u = 2:4
    ok = ok && (bitxor(got(1), got(u)) == bitxor(block(1), block(u)));
end
results(end+1,:) = {ok, 'T3: XOR preservation (Eq. 8)'};

%% ── T4: CRT-P pixel reconstruction ──────────────────────────────────────────
fprintf('--- T4: CRT-P reconstruction ---\n');
rng(7);
ok = true;
combos3 = nchoosek(1:4, 3);
for t = 1:300
    m  = uint64(randi([0, 255]));
    r1 = uint64(randi([0, 255]));
    r2 = uint64(randi([0, 65535]));
    sh = zeros(1, 4, 'uint64');
    for i = 1:4
        Y = bitxor(bitxor(m, r1), gf2_mul(r2, q0));
        sh(i) = gf2_mod(Y, q_list(i));
    end
    for c = 1:size(combos3, 1)
        idx = combos3(c, :);
        Q = uint64(1);
        for i = idx
            Q = gf2_mul(Q, q_list(i));
        end
        Yr = uint64(0);
        for i = idx
            [Qi, ~] = gf2_divmod(Q, q_list(i));
            inv = gf2_inv(gf2_mod(Qi, q_list(i)), q_list(i));
            Yr = bitxor(Yr, gf2_mod(gf2_mul(gf2_mul(sh(i), Qi), inv), Q));
        end
        mr = bitxor(gf2_mod(Yr, q0), r1);
        if mr ~= m
            ok = false;
        end
    end
end
results(end+1,:) = {ok, 'T4: 300 pixels x all C(4,3) combos reconstruct'};

%% ── T5: scramble inverse ────────────────────────────────────────────────────
fprintf('--- T5: scramble/unscramble ---\n');
rng(5);
img5 = uint8(randi([0, 255], 32, 48));
scr5 = scramble_image(img5, [11, 22]);
results(end+1,:) = {isequal(unscramble_image(scr5, [11, 22]), img5), ...
    'T5: unscramble(scramble(img)) == img'};
results(end+1,:) = {~isequal(scr5, img5), 'T5: scrambling changes the image'};

%% ── T6: AC embed/extract round-trip ─────────────────────────────────────────
fprintf('--- T6: AC embed/extract round-trip ---\n');
rng(6);
img6 = make_test_gray_image('smooth', 32, 32);
scr6 = scramble_image(img6, [1, 2]);
sh6  = crtp_share(scr6, 4, 3, 3, 4, q_list, q0);
ok_share = true; ok_sec = true;
for i = 1:4
    [~, ~, st] = ac_embed(sh6{i}, [], 50 + i);
    sec = randi([0, 1], 1, st.pureEC);          % fill to FULL pure capacity
    [mk, sd] = ac_embed(sh6{i}, sec, 50 + i);
    [rest, sec_out] = ac_extract(mk, sd, 50 + i);
    ok_share = ok_share && isequal(rest, sh6{i});
    ok_sec   = ok_sec   && isequal(sec_out, sec);
end
results(end+1,:) = {ok_share, 'T6: all 4 shares restored losslessly'};
results(end+1,:) = {ok_sec,   'T6: secrets recovered bit-exact at full pure EC'};

%% ── T7: full pipeline, all C(4,3) combos ────────────────────────────────────
fprintf('--- T7: full pipeline ---\n');
types = {'smooth', 'mixed'};
for tt = 1:numel(types)
    rng(70 + tt);
    img7 = make_test_gray_image(types{tt}, 32, 32);
    scr7 = scramble_image(img7, [9, 8]);
    sh7  = crtp_share(scr7, 4, 3, 7, 6, q_list, q0);
    mk7 = cell(1, 4); sd7 = cell(1, 4); sec7 = cell(1, 4);
    for i = 1:4
        [~, ~, st] = ac_embed(sh7{i}, [], 90 + i);
        sec7{i} = randi([0, 1], 1, max(0, st.pureEC - 5));
        [mk7{i}, sd7{i}] = ac_embed(sh7{i}, sec7{i}, 90 + i);
    end
    s_total = 32 * 32 / 4;
    ok = true;
    for c = 1:size(combos3, 1)
        idx = combos3(c, :);
        rest = cell(1, 3);
        for j = 1:3
            [rest{j}, so] = ac_extract(mk7{idx(j)}, sd7{idx(j)}, 90 + idx(j));
            ok = ok && isequal(so, sec7{idx(j)});
        end
        rec = unscramble_image( ...
            crtp_reconstruct(rest, q_list(idx), q0, 7, s_total), [9, 8]);
        ok = ok && isequal(rec, img7);
    end
    results(end+1,:) = {ok, ...
        sprintf('T7: %s image — all C(4,3) combos lossless + secrets exact', types{tt})};
end

%% ── T8: complex (B2-heavy) image ────────────────────────────────────────────
fprintf('--- T8: complex image ---\n');
rng(88);
img8 = make_test_gray_image('complex', 32, 32);
scr8 = scramble_image(img8, [3, 4]);
sh8  = crtp_share(scr8, 4, 3, 5, 6, q_list, q0);
[~, ~, st8] = ac_embed(sh8{1}, [], 99);
sec8 = randi([0, 1], 1, max(0, st8.pureEC));
[mk8, sd8] = ac_embed(sh8{1}, sec8, 99);
[rest8, sec8out] = ac_extract(mk8, sd8, 99);
results(end+1,:) = {isequal(rest8, sh8{1}) && isequal(sec8out, sec8), ...
    sprintf('T8: complex image restores (pure EC=%d, %d B2 blocks)', ...
            st8.pureEC, st8.n_B2)};

%% ── T9: capacity overflow raises an error ───────────────────────────────────
fprintf('--- T9: capacity overflow ---\n');
overflow_caught = false;
try
    ac_embed(sh8{1}, ones(1, 32*32*8), 99);   % absurdly large payload
catch
    overflow_caught = true;
end
results(end+1,:) = {overflow_caught, 'T9: oversized payload raises error'};

%% ── Summary ─────────────────────────────────────────────────────────────────
n_pass = 0; n_fail = 0;
fprintf('\n');
for r = 1:size(results, 1)
    if results{r, 1}
        fprintf('  [PASS] %s\n', results{r, 2});
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] %s\n', results{r, 2});
        n_fail = n_fail + 1;
    end
end
fprintf('\n========================================\n');
fprintf('  Results: %d passed, %d failed\n', n_pass, n_fail);
fprintf('========================================\n');
if n_fail == 0
    fprintf('  ALL TESTS PASSED\n');
else
    fprintf('  SOME TESTS FAILED — check output above\n');
end
