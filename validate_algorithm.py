"""Algorithm validation for the CRTP-SS + adaptive-coding RDHEI scheme
(Fang et al., IEEE IoTJ 2025) before committing to the MATLAB port.

Validates:
  V1: GF(2)[x] arithmetic (mul, mod, extended-Euclid inverse)
  V2: Paper's worked example (Sec. IV-A, Eqs. 9-12): block [192 191; 195 180]
      with R1=10000000, R2=...0001, q0=111110101, q1=100011011
      must encrypt to share1 block [174 209; 173 218]
  V3: XOR-preservation property (Eq. 8): me_1 ^ me_u == m_1 ^ m_u
  V4: CRT-P reconstruction (Eq. 14) from any k=3 of n=4 shares
  V5: Adaptive-coding embed -> extract round-trip on a synthetic image
  V6: Full pipeline: scramble -> share -> embed -> extract -> reconstruct ->
      unscramble == original, and secret recovered bit-exact
"""
import random

# ---------------- GF(2)[x] arithmetic ----------------
def gf2_deg(a):
    return a.bit_length() - 1

def gf2_mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r

def gf2_divmod(a, b):
    q = 0
    db = gf2_deg(b)
    while a and gf2_deg(a) >= db:
        sh = gf2_deg(a) - db
        q ^= 1 << sh
        a ^= b << sh
    return q, a

def gf2_mod(a, b):
    return gf2_divmod(a, b)[1]

def gf2_inv(a, m):
    # extended Euclid: find t with a*t = 1 mod m
    r0, r1 = m, gf2_mod(a, m)
    t0, t1 = 0, 1
    while r1:
        q, r = gf2_divmod(r0, r1)
        r0, r1 = r1, r
        t0, t1 = t1, t0 ^ gf2_mul(q, t1)
    assert r0 == 1, "not invertible"
    return t0

# V1 sanity: inverse property in GF(2^8) with AES poly
q_aes = 0b100011011
for a in [1, 2, 3, 0x53, 0xCA, 0xFF]:
    inv = gf2_inv(a, q_aes)
    assert gf2_mod(gf2_mul(a, inv), q_aes) == 1
print("V1 PASS: GF(2)[x] mul/mod/inv")

# ---------------- V2: paper worked example ----------------
q0 = 0b111110101
q1 = 0b100011011
q2 = 0b100011101
q3 = 0b100101011
q4 = 0b100101101
R1 = 0b10000000
R2 = 0b0000000000000001

def crtp_encrypt_pixel(m, R1z, R2z, q0, qi):
    Y = m ^ R1z ^ gf2_mul(R2z, q0)
    return gf2_mod(Y, qi)

block = [192, 191, 195, 180]
share1 = [crtp_encrypt_pixel(m, R1, R2, q0, q1) for m in block]
assert share1 == [174, 209, 173, 218], f"expected [174,209,173,218], got {share1}"
print(f"V2 PASS: worked example block {block} -> share1 {share1}")

# ---------------- V3: XOR preservation ----------------
for u in range(1, 4):
    assert (share1[0] ^ share1[u]) == (block[0] ^ block[u])
print("V3 PASS: XOR preservation (Eq. 8)")

# ---------------- V4: CRT-P reconstruction ----------------
def crtp_reconstruct_pixel(shares, qs, q0, R1z):
    # shares: list of k share pixels; qs: their moduli
    Q = 1
    for q in qs:
        Q = gf2_mul(Q, q)
    Y = 0
    for me, q in zip(shares, qs):
        Qi, rem = gf2_divmod(Q, q)
        assert rem == 0
        inv = gf2_inv(gf2_mod(Qi, q), q)
        Y ^= gf2_mod(gf2_mul(gf2_mul(me, Qi), inv), Q)
    return gf2_mod(Y, q0) ^ R1z

qs_all = [q1, q2, q3, q4]
random.seed(7)
ok = True
for trial in range(2000):
    m = random.randrange(256)
    r1 = random.randrange(256)
    r2 = random.randrange(65536)
    shares = [crtp_encrypt_pixel(m, r1, r2, q0, q) for q in qs_all]
    # every k=3 combination of the 4 shares
    import itertools
    for idx in itertools.combinations(range(4), 3):
        rec = crtp_reconstruct_pixel([shares[i] for i in idx],
                                     [qs_all[i] for i in idx], q0, r1)
        if rec != m:
            ok = False
            print(f"FAIL m={m} r1={r1} r2={r2} idx={idx} rec={rec}")
assert ok
print("V4 PASS: CRT-P reconstruction, 2000 pixels x all C(4,3) combos")

# ---------------- V5/V6: adaptive coding + full pipeline ----------------
import heapq

def classify_lambda(Dmax):
    if Dmax == 0:
        return 0
    return max(1, Dmax.bit_length())  # ceil(log2(Dmax+1))

def huffman_codes(weights):
    """weights: dict sym->w (w>0). Returns dict sym->code string."""
    if len(weights) == 1:
        return {next(iter(weights)): "0"}
    h = [(w, i, (s,)) for i, (s, w) in enumerate(sorted(weights.items()))]
    heapq.heapify(h)
    codes = {s: "" for s in weights}
    cnt = len(h)
    while len(h) > 1:
        w1, _, g1 = heapq.heappop(h)
        w2, _, g2 = heapq.heappop(h)
        for s in g1:
            codes[s] = "0" + codes[s]
        for s in g2:
            codes[s] = "1" + codes[s]
        cnt += 1
        heapq.heappush(h, (w1 + w2, cnt, g1 + g2))
    return codes

def build_label_codes(lam_hist):
    """Huffman over lambda histogram with merge rule: while max code len > 5,
    merge min-weight lambda into the adjacent lambda (toward larger lambda
    when possible) so its blocks are treated more conservatively."""
    merge_map = {l: l for l in range(7)}
    while True:
        eff = {}
        for l, w in lam_hist.items():
            if w > 0:
                eff[merge_map[l]] = eff.get(merge_map[l], 0) + w
        codes = huffman_codes(eff)
        if not codes or max(len(c) for c in codes.values()) <= 5:
            return codes, merge_map
        # merge the min-weight effective symbol into next higher (or lower for 6)
        lmin = min(eff, key=lambda s: (eff[s], s))
        target = lmin + 1 if lmin < 6 else 5
        # NOTE: merging lambda L into L-1 is UNSAFE (d won't fit);
        # for lmin==6 we instead merge 5 into 6.
        if lmin == 6:
            src, dst = 5, 6
        else:
            src, dst = lmin, lmin + 1
        for l in range(7):
            if merge_map[l] == src:
                merge_map[l] = dst

def ac_embed(share, M, N, payload_secret):
    """share: dict (r,c)->pixel of an encrypted share, 2x2 blocks.
    Returns marked share, side_info. Payload = labels + B2 MSBs + secret."""
    blocks = []
    for br in range(0, M, 2):
        for bc in range(0, N, 2):
            p = [share[(br, bc)], share[(br, bc+1)], share[(br+1, bc)], share[(br+1, bc+1)]]
            D = [p[0] ^ p[u] for u in (1, 2, 3)]
            Dmax = max(D)
            if Dmax <= 63:
                lam = classify_lambda(Dmax)
                blocks.append(("B1", (br, bc), p, D, lam))
            else:
                blocks.append(("B2", (br, bc), p, D, None))

    lam_hist = {}
    for b in blocks:
        if b[0] == "B1":
            lam_hist[b[4]] = lam_hist.get(b[4], 0) + 1
    codes, merge_map = build_label_codes(lam_hist) if lam_hist else ({}, {l: l for l in range(7)})

    b1 = [b for b in blocks if b[0] == "B1"]
    label_bits = "".join(codes[merge_map[b[4]]] for b in b1[1:])
    msb_bits = "".join(str((b[2][0] >> 7) & 1) for b in blocks if b[0] == "B2")
    stream = label_bits + msb_bits + payload_secret

    caps = [3 * (8 - merge_map[b[4]]) - 1 for b in b1]
    assert len(stream) <= sum(caps), "payload exceeds capacity"
    stream += "0" * (sum(caps) - len(stream))  # pad

    marked = dict(share)
    pos = 0
    for (typ, (br, bc), p, D, lam) in blocks:
        if typ == "B2":
            marked[(br, bc)] = (p[0] & 0x7F) | 0x80  # type bit 1
            continue
        lam_e = merge_map[lam]
        e1 = lam_e
        e2 = 8 - e1
        marked[(br, bc)] = p[0] & 0x7F  # type bit 0
        coords = [(br, bc+1), (br+1, bc), (br+1, bc+1)]
        for ui, (r, c) in enumerate(coords):
            if ui == 0:
                head = str((p[0] >> 7) & 1) + stream[pos:pos + e2 - 1]
                pos += e2 - 1
            else:
                head = stream[pos:pos + e2]
                pos += e2
            dbits = format(D[ui], f"0{e1}b") if e1 else ""
            marked[(r, c)] = int(head + dbits, 2)
    side = {
        "lambda_first": merge_map[b1[0][4]] if b1 else None,
        "codes": codes,
        "n_secret": len(payload_secret),
    }
    return marked, side

def ac_extract(marked, M, N, side):
    codes = side["codes"]
    dec = {v: k for k, v in codes.items()}
    restored = dict(marked)
    buf = []          # extracted payload bits
    label_ptr = 0
    b1_seen = 0
    b2_blocks = []
    for br in range(0, M, 2):
        for bc in range(0, N, 2):
            ref = marked[(br, bc)]
            if ref >> 7 == 1:   # B2
                b2_blocks.append((br, bc))
                continue
            b1_seen += 1
            if b1_seen == 1:
                lam = side["lambda_first"]
            else:
                code = ""
                while code not in dec:
                    assert label_ptr < len(buf), "label bits not yet available!"
                    code += buf[label_ptr]
                    label_ptr += 1
                lam = dec[code]
            e1, e2 = lam, 8 - lam
            coords = [(br, bc+1), (br+1, bc), (br+1, bc+1)]
            Ds = []
            ref_msb = None
            for ui, (r, c) in enumerate(coords):
                v = format(marked[(r, c)], "08b")
                if ui == 0:
                    ref_msb = v[0]
                    buf.extend(v[1:e2])
                else:
                    buf.extend(v[:e2])
                Ds.append(int(v[e2:], 2) if e1 else 0)
            ref_orig = int(ref_msb) * 128 + (ref & 0x7F)
            restored[(br, bc)] = ref_orig
            for ui, (r, c) in enumerate(coords):
                restored[(r, c)] = ref_orig ^ Ds[ui]
    # consume stream: labels already consumed up to label_ptr
    rest = buf[label_ptr:]
    for i, (br, bc) in enumerate(b2_blocks):
        msb = int(rest[i])
        restored[(br, bc)] = msb * 128 + (marked[(br, bc)] & 0x7F)
    secret = "".join(rest[len(b2_blocks):len(b2_blocks) + side["n_secret"]])
    return restored, secret

# Build a synthetic smooth-ish grayscale image (mix of smooth + noisy regions)
random.seed(123)
M, N = 16, 16
img = {}
for r in range(M):
    for c in range(N):
        if c < N // 2:
            img[(r, c)] = min(255, 100 + r * 3 + (c % 4))       # smooth
        else:
            img[(r, c)] = random.randrange(256)                   # complex

# scramble blocks + pixels (Key0)
rs = random.Random(999)
bpos = [(br, bc) for br in range(0, M, 2) for bc in range(0, N, 2)]
perm = bpos[:]
rs.shuffle(perm)
inblock = {b: rs.sample(range(4), 4) for b in bpos}
off = [(0, 0), (0, 1), (1, 0), (1, 1)]
scr = {}
for (src, dst) in zip(bpos, perm):
    pp = inblock[src]
    for u in range(4):
        s_off = off[pp[u]]
        d_off = off[u]
        scr[(dst[0] + d_off[0], dst[1] + d_off[1])] = img[(src[0] + s_off[0], src[1] + s_off[1])]

# CRTP-SS share (Key1, Key2) — per-block R1, R2
rk = random.Random(555)
R1z = {b: rk.randrange(256) for b in bpos}
R2z = {b: rk.randrange(65536) for b in bpos}
n_sh, k_sh = 4, 3
shares = [dict() for _ in range(n_sh)]
for b in bpos:
    for u in range(4):
        pos = (b[0] + off[u][0], b[1] + off[u][1])
        for i, q in enumerate(qs_all):
            shares[i][pos] = crtp_encrypt_pixel(scr[pos], R1z[b], R2z[b], q0, q)

# embed different secret in each share
secrets = ["".join(rk.choice("01") for _ in range(40)) for _ in range(n_sh)]
marked, sides = [], []
for i in range(n_sh):
    mk, sd = ac_embed(shares[i], M, N, secrets[i])
    marked.append(mk)
    sides.append(sd)

# extract + restore from each share, verify secrets and share restoration
for i in range(n_sh):
    rest, sec = ac_extract(marked[i], M, N, sides[i])
    assert sec == secrets[i], f"secret mismatch share {i}"
    assert rest == shares[i], f"share {i} not restored losslessly"
print("V5 PASS: AC embed/extract round-trip, all 4 shares restored + secrets exact")

# reconstruct from k=3 restored shares, unscramble, compare
import itertools
for idx in itertools.combinations(range(4), 3):
    rec_scr = {}
    for b in bpos:
        for u in range(4):
            pos = (b[0] + off[u][0], b[1] + off[u][1])
            rec_scr[pos] = crtp_reconstruct_pixel(
                [shares[i][pos] for i in idx], [qs_all[i] for i in idx], q0, R1z[b])
    # unscramble
    rec = {}
    for (src, dst) in zip(bpos, perm):
        pp = inblock[src]
        for u in range(4):
            s_off = off[pp[u]]
            d_off = off[u]
            rec[(src[0] + s_off[0], src[1] + s_off[1])] = rec_scr[(dst[0] + d_off[0], dst[1] + d_off[1])]
    assert rec == img, f"recovery failed for combo {idx}"
print("V6 PASS: full pipeline lossless for all C(4,3) share combos")
print("\nALL VALIDATIONS PASSED — algorithm design is correct.")
