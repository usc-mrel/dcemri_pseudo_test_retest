function kern = VSkern(type,fw,fwhm)
%   Generates a view sharing kernel

if ~mod(fw,2)
    fw = fw+1;
end

%   Compute window
switch lower(type)
    case {'rect','rectwin'}
        kern = ones([fw,1]);
    
    case {'gauss','gausswin'}
        kern = exp(-((1:fw)-ceil(fw/2)).^2/(2*(fwhm/(2*sqrt(2*log(2)))).^2));
    
    otherwise
        error('Unkown view-sharing kernel');
end

%   Zero the center value and normalize
kern(ceil(fw/2)) = 0;
kern = kern ./ sum(kern(:));

end

