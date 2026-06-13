function r = gf2_mul(a, b)
% GF2_MUL  Carry-less multiplication of GF(2) polynomials (uint64 bitmasks).
%          deg(a) + deg(b) must be <= 63 (caller's responsibility).
a = uint64(a);
b = uint64(b);
r = uint64(0);
while b > 0
    if bitand(b, uint64(1)) > 0
        r = bitxor(r, a);
    end
    a = bitshift(a, 1);
    b = bitshift(b, -1);
end
end
