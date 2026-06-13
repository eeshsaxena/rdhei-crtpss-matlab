function t = gf2_inv(a, m)
% GF2_INV  Modular inverse of polynomial a modulo m over GF(2) via
%          the extended Euclidean algorithm: a * t == 1 (mod m).
r0 = uint64(m);
r1 = gf2_mod(a, m);
t0 = uint64(0);
t1 = uint64(1);
while r1 ~= 0
    [q, r] = gf2_divmod(r0, r1);
    r0 = r1;  r1 = r;
    tn = bitxor(t0, gf2_mul(q, t1));
    t0 = t1;  t1 = tn;
end
if r0 ~= 1
    error('gf2_inv: polynomial is not invertible modulo m (gcd != 1)');
end
t = t0;
end
