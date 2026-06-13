function d = gf2_deg(a)
% GF2_DEG  Degree of a GF(2) polynomial stored as a uint64 bitmask.
%          gf2_deg(0) = -1 by convention.
a = uint64(a);
d = -1;
while a > 0
    a = bitshift(a, -1);
    d = d + 1;
end
end
