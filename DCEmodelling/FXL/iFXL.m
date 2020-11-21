function C = iFXL(R1, opt)
% inverse Fast-Exchange-Limit
%
% Yannick Bliesener 2019

% compatibility with older Matlab versions
if numel(opt.R1) == 1
    C = R1 - opt.R1;
else
    rvec = [ones(1,ndims(R1)-1), size(R1,ndims(R1))];
    C = R1 - repmat(opt.R1, rvec);
end
C = C / opt.TK.r1;


end

