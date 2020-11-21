function varargout = iSPGR(img, opt)
% inverse SPGR sequence model
%
%   examples:
%       1) R1 = iSPGR(img, opt)
%       2) [R1, deltaB] = iSPGR(img, opt)
%       3) [R1, R2, deltaB] = iSPGR(img, opt)
%
%   Yannick 2019

alpha = opt.FA * pi / 180;

%   B1
if isfield(opt,'B1') && ~isempty(opt.B1)
    alpha = alpha .* opt.B1;
end

switch nargout
    
    case 1
        img = real(img);
        img( img(:) < 0 ) = 0;
        
        varargout{1} = sig2R1(img, alpha, opt);
        
    case 2
        varargout{1} = sig2R1(abs(img), alpha, opt);
        varargout{2} = sig2deltaB(img, opt);
    case 3
        error('TODO')
        % this one would require knowledge of R2/R2star
    otherwise
        error('Unkown parameterization')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Utilities %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function R1 = sig2R1(img, alpha, opt)
        
        rvec = [ones(1,ndims(img)-1), size(img,ndims(img))];
            % compatibility with older Matlab versions
        
        A = img ./ repmat(( (opt.M0 + eps) .* sin( alpha) ), rvec);
        
        A( repmat((opt.M0 < eps),rvec) ) = 0;
            % correct for numerical instability if there is no tissue at a
            % certain location

        cosa = cos( alpha );
        if numel(cosa) > 1
            % compatibility with older Matlab versions
            E = ( A - 1 ) ./ ( A.* repmat(cosa, rvec) - 1 );
        else
            E = ( A - 1 ) ./ ( A.* cosa - 1 );
        end

        indError = (E <= 0.0) | (E > 1.0);

        if any( indError )
            % this probably happens when the tissue is non-enhancing, and there is
            % noise so the difference between time series and baseline is negative
            % => make sure that tissue is set to zeros concentration!

            warning('Error in sig2conc: Cannot convert signal to concentration!')
            fprintf('Correcting %i violation(s) in E1!\n', sum(indError(:)) );
        end
        
        R1 = log( E ) / ( -opt.TR );
        R1 = real(R1);
        
        R1( E(:) <= 0.0 ) = Inf;
        R1( E(:) > 1.0 ) = nan;
    end

    function deltaB = sig2deltaB(img, opt)
        deltaB = angle(img);
        
        s = size(deltaB);
        deltaB = reshape(deltaB, [], s(end));
        for i=1:size(deltaB,1)
            deltaB(i,:) = unwrap(deltaB(i,:));
        end
        deltaB = reshape(deltaB, s);
        
        deltaB = deltaB / (-opt.omega0*opt.TE);
    end

end
