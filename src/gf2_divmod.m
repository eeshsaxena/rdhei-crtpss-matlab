function [q, r] = gf2_divmod(a, b)
% GF2_DIVMOD  Polynomial long division over GF(2): a = q*b + r, deg(r) < deg(b).
a = uint64(a);
b = uint64(b);
if b == 0
    error('gf2_divmod: division by zero polynomial');
end
q = uint64(0);
db = gf2_deg(b);
while a ~= 0 && gf2_deg(a) >= db
    sh = gf2_deg(a) - db;
    q = bitxor(q, bitshift(uint64(1), sh));
    a = bitxor(a, bitshift(b, sh));
end
r = a;
end
