function [R1, R2] = fFXL(C, opt)
% forward Fast-Exchange-Limit
%
% Yannick Bliesener 2019

rvec = [ones(1,ndims(C)-1), size(C,ndims(C))];
    % compatibility with older Matlab versions

if numel(opt.R1) == 1
    R1 = opt.R1 + opt.TK.r1 * C;
else
    R1 = repmat(opt.R1, rvec) + opt.TK.r1 * C;
end

if nargout > 1
    if numel(opt.R2) == 1
        R2 = opt.R2 + opt.TK.r2 * C;
    else
        R2 = repmat(opt.R2, rvec) + opt.TK.r2 * C;
    end
end

end

