function r = gf2_mod(a, b)
% GF2_MOD  Remainder of GF(2) polynomial division (uint64 bitmasks).
[~, r] = gf2_divmod(a, b);
end
