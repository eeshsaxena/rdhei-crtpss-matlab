function [q0, q_list, all30] = irreducible_polys(n)
% IRREDUCIBLE_POLYS  Irreducible degree-8 polynomials over GF(2) used as
%                    CRTP-SS moduli (Sec. III-B / Sec. V-A of the paper).
%
%   [q0, q_list, all30] = irreducible_polys(n)
%
%   q0     : auxiliary modulus q0(x) = 111110101b   (paper Sec. V-A)
%   q_list : 1 x n vector of share moduli q1..qn (paper's values for n<=4,
%            then further irreducible polynomials from the canonical list)
%   all30  : all 30 irreducible degree-8 polynomials over GF(2)
%            (key-space basis, Eq. 16)
%
%   Distinct irreducible polynomials are automatically pairwise coprime,
%   satisfying the CRTP-SS constraints gcd(qi,qj)=1 and gcd(qi,q0)=1.

% Paper's experimental parameters (Sec. V-A):
%   q0 = 111110101, q1 = 100011011, q2 = 100011101,
%   q3 = 100101011, q4 = 100101101
q0 = uint64(bin2dec('111110101'));      % 0x1F5

paper_q = uint64([ ...
    bin2dec('100011011'), ...           % 0x11B (AES polynomial)
    bin2dec('100011101'), ...           % 0x11D
    bin2dec('100101011'), ...           % 0x12B
    bin2dec('100101101')]);             % 0x12D

% The 30 irreducible degree-8 polynomials over GF(2), as 9-bit masks.
all30 = uint64(hex2dec({ ...
    '11B','11D','12B','12D','139','13F','14D','15F','163','165', ...
    '169','171','177','17B','187','18B','18D','19F','1A3','1A9', ...
    '1B1','1BD','1C3','1CF','1D7','1DD','1E7','1F3','1F5','1F9'}))';

if n > 29
    error('At most 29 share moduli available (one of the 30 is used as q0).');
end

if n <= numel(paper_q)
    q_list = paper_q(1:n);
else
    % Extend with further irreducible polys, skipping q0 and paper's set
    pool = setdiff(all30, [q0, paper_q], 'stable');
    q_list = [paper_q, pool(1:(n - numel(paper_q)))];
end
end
