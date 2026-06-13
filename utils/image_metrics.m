function out = image_metrics(mode, A, B)
% IMAGE_METRICS  Security/quality metrics used in Sec. V of the paper.
%
%   H    = image_metrics('entropy', A)       % information entropy (Eq. 17)
%   v    = image_metrics('npcr',    A, B)    % NPCR in percent (ideal 99.6094)
%   v    = image_metrics('uaci',    A, B)    % UACI in percent (ideal 33.4635)
%   v    = image_metrics('psnr',    A, B)    % PSNR in dB (Inf if identical)

switch lower(mode)
    case 'entropy'
        h = histcounts(double(A(:)), -0.5:1:255.5);
        p = h / sum(h);
        p = p(p > 0);
        out = -sum(p .* log2(p));

    case 'npcr'
        out = mean(A(:) ~= B(:)) * 100;

    case 'uaci'
        out = mean(abs(double(A(:)) - double(B(:))) / 255) * 100;

    case 'psnr'
        mse = mean((double(A(:)) - double(B(:))).^2);
        if mse == 0
            out = Inf;
        else
            out = 10 * log10(255^2 / mse);
        end

    otherwise
        error('Unknown metric: %s', mode);
end
end
