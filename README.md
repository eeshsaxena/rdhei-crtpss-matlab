# RDHEI with CRTP-SS and Adaptive Coding — MATLAB Implementation

MATLAB implementation of:

> G. Fang, F. Wang, C. Zhao, C. Qin, C.-C. Chang, C.-C. Chang,
> **"Reversible Data Hiding With Secret Encrypted Image Sharing and Adaptive Coding"**
> *IEEE Internet of Things Journal*, vol. 12, no. 13, pp. 23933–23945, July 2025.
> DOI: 10.1109/JIOT.2025.3555380

---

## Overview

A multi-party RDHEI scheme with two key components:

1. **CRTP-SS encryption** — secret sharing based on the Chinese Remainder
   Theorem for **polynomials** over GF(2⁸), using irreducible degree-8
   polynomials as moduli. Unlike integer-CRT schemes, encrypted pixels are
   uniformly distributed over the full 0–255 range, and the in-block XOR
   structure is preserved (Eq. 8), vacating room for embedding with **no
   preprocessing**.
2. **Adaptive-coding (AC) embedding** — XOR-preservation (XORP) compresses
   each 2×2 block against its reference pixel; per-block thresholds λ are
   labeled with **Huffman codes** (≤ 5 bits, guaranteed decodable on the fly).

The original image is recovered **losslessly** from any *k* of *n* marked
shares; secrets embed/extract bit-exact per data hider.

```
Content owner : scramble (Key0) ──> CRTP-SS (k,n) ──> n encrypted shares
Data hider i  : AC embed (Key3_i) ──> marked share i
Receiver      : any k shares ──> AC extract + restore ──> CRT-P ──> unscramble
```

## Repo Structure

| File | Implements |
|------|------------|
| `src/gf2_{deg,mul,divmod,mod,inv}.m` | GF(2)[x] polynomial arithmetic (uint64 bitmasks) |
| `src/irreducible_polys.m` | Paper's moduli q₀…q₄ + all 30 irreducible degree-8 polys |
| `src/scramble_image.m` / `unscramble_image.m` | Key₀ block + pixel scrambling (Sec. IV-A) |
| `src/crtp_share.m` | CRTP-SS encryption, Eqs. 5–6 |
| `src/crtp_reconstruct.m` | CRT-P reconstruction, Eq. 14 |
| `src/huffman_label.m` | λ label codes + threshold-merge rule (Sec. IV-B) |
| `src/ac_embed.m` | Adaptive-coding embedding (Fig. 4) |
| `src/ac_extract.m` | Extraction + share restoration (Fig. 5) |
| `utils/image_metrics.m` | Entropy (Eq. 17), NPCR, UACI, PSNR |
| `main_demo.m` | Full (3,4)-threshold pipeline demo |
| `run_all_tests.m` | 9 test groups incl. the paper's worked example |
| `run_experiments.m` | Tables II/VI/VII + Fig. 10 equivalents |
| `validate_algorithm.py` | Python pre-validation of the algorithm design |

## Quick Start

```matlab
addpath(genpath('src')); addpath(genpath('utils'));

k = 3; n = 4;
[q0, q_list] = irreducible_polys(n);
img = imread('house.png');                 % 8-bit grayscale, even dims

% Content owner
scr    = scramble_image(img, [11 22]);     % Key0
shares = crtp_share(scr, n, k, 33, 44, q_list, q0);   % Key1, Key2

% Data hider i (independent)
[~, ~, st] = ac_embed(shares{1}, [], 101); % query capacity
secret = randi([0 1], 1, st.pureEC);
[marked, side] = ac_embed(shares{1}, secret, 101);    % Key3

% Receiver (any k shares)
[restored, secret_out] = ac_extract(marked, side, 101);
% ... collect k restored shares, then:
s_total = numel(img) / 4;
rec = unscramble_image( ...
        crtp_reconstruct({restored, r2, r3}, q_list([1 2 3]), q0, 33, s_total), ...
        [11 22]);
isequal(rec, img)        % 1 — lossless
```

Run `main_demo`, `run_all_tests`, or `run_experiments` directly.

## Verified Against the Paper

`run_all_tests.m` T2 checks the paper's worked example (Sec. IV-A):
block `[192 191; 195 180]` with R1=`10000000`, R2=`…0001`,
q₀=`111110101`, q₁=`100011011` encrypts to share-1 block
`[174 209; 173 218]` — bit-exact. The same example and the full pipeline
were independently validated in `validate_algorithm.py` (all 6 checks pass,
including 2000 random pixels × all C(4,3) reconstruction combos).

## Algorithm Notes

- **λ classification** (Table I): λ = 0 if Dmax = 0, else ⌈log₂(Dmax+1)⌉;
  B1 iff Dmax ≤ 63. Capacity per B1 block = 3·(8−λ) − 1 bits.
- **Label decodability**: max Huffman code length is forced ≤ 5 bits
  (= minimum B1 capacity), so the receiver always holds block *j*'s label
  before reaching block *j*. When the tree depth would be 6 (7 active
  thresholds), the min-weight threshold is merged into the next one.
- **R2 width** = min(16, 8(k−1)) bits so deg Y < 8k, guaranteeing a unique
  CRT solution (paper's (3,4) setting → 16 bits).
- **Side info** (`side_info` struct): Huffman table, first-block λ, secret
  length — the paper's "informed" auxiliary information.

## Requirements

MATLAB R2019b+ (or recent Octave). No toolboxes required.

## Companion Repo

[mrdhcbi-matlab](../mrdhcbi-matlab) — Chen et al., *"Multi-Party Reversible
Data Hiding in Ciphertext Binary Images Based on Visual Cryptography"*,
IEEE SPL 2025.

## Reference

```bibtex
@article{fang2025rdhei,
  author  = {Fang, Guangtian and Wang, Feng and Zhao, Chenbin and Qin, Chuan
             and Chang, Ching-Chun and Chang, Chin-Chen},
  title   = {Reversible Data Hiding With Secret Encrypted Image Sharing
             and Adaptive Coding},
  journal = {IEEE Internet of Things Journal},
  year    = {2025},
  volume  = {12},
  number  = {13},
  pages   = {23933--23945},
  doi     = {10.1109/JIOT.2025.3555380}
}
```
