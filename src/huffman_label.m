function [codes, merge_map] = huffman_label(lam_hist)
% HUFFMAN_LABEL  Build Huffman label codes for thresholds lambda = 0..6
%                (Sec. IV-B2, Fig. 3, Table I) with the threshold-merging
%                rule of Sec. IV-B3.
%
%   [codes, merge_map] = huffman_label(lam_hist)
%
%   lam_hist  : 1x7 vector, lam_hist(l+1) = number of B1 blocks with
%               threshold lambda = l
%   codes     : 1x7 cell array of char codes; codes{l+1} = label code for
%               effective lambda l ('' if l never used as effective value)
%   merge_map : 1x7 vector, merge_map(l+1) = effective lambda for blocks
%               classified at lambda l (after merging)
%
%   Merge rule: the maximum label code length must be <= 5 bits so the
%   receiver always has the next block's label available (min B1 capacity
%   is 3*(8-6)-1 = 5 bits). While the Huffman tree is deeper than 5, the
%   minimum-weight threshold is merged into the NEXT (higher) threshold —
%   safe because a larger lambda stores d values in more bits. If the
%   minimum is lambda=6, threshold 5 is merged into 6 instead.

merge_map = 0:6;

while true
    % Effective histogram after merging
    eff = zeros(1, 7);
    for l = 0:6
        if lam_hist(l+1) > 0
            e = merge_map(l+1);
            eff(e+1) = eff(e+1) + lam_hist(l+1);
        end
    end
    syms = find(eff > 0) - 1;      % effective lambda values present
    w    = eff(syms + 1);

    codes = repmat({''}, 1, 7);
    if isempty(syms)
        return;                    % no embeddable blocks at all
    end
    if numel(syms) == 1
        codes{syms+1} = '0';
        return;
    end

    % --- standard Huffman over (syms, w), deterministic tie-breaking ---
    % Each tree node: weight + member list; groups stored in cell array.
    grp_w = num2cell(w);
    grp_m = num2cell(syms);
    code_acc = repmat({''}, 1, 7);
    while numel(grp_w) > 1
        ws = cell2mat(grp_w);
        [~, order] = sortrows([ws(:), cellfun(@min, grp_m(:))]);
        i1 = order(1); i2 = order(2);
        for mlist = grp_m{i1}
            code_acc{mlist+1} = ['0', code_acc{mlist+1}];
        end
        for mlist = grp_m{i2}
            code_acc{mlist+1} = ['1', code_acc{mlist+1}];
        end
        merged_w = grp_w{i1} + grp_w{i2};
        merged_m = [grp_m{i1}, grp_m{i2}];
        keep = setdiff(1:numel(grp_w), [i1, i2]);
        grp_w = [grp_w(keep), {merged_w}];
        grp_m = [grp_m(keep), {merged_m}];
    end
    for l = syms
        codes{l+1} = code_acc{l+1};
    end

    maxlen = max(cellfun(@length, codes(syms+1)));
    if maxlen <= 5
        return;
    end

    % --- merge rule: min-weight effective threshold -> next threshold ---
    [~, imin] = min(w);
    lmin = syms(imin);
    if lmin == 6
        src = 5; dst = 6;          % can't go above 6; fold 5 into 6
    else
        src = lmin; dst = lmin + 1;
    end
    for l = 0:6
        if merge_map(l+1) == src
            merge_map(l+1) = dst;
        end
    end
end
end
