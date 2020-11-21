function varargout = vecsplit(x)
%function varargout = vecsplit(x)
%   usage:
%       dims = [Nx, Ny, Nz, Nt];
%       [Nx,Ny,Nz,Nt] = vecsplit(dims);
%   Yannick 2016

nout = max(nargout,1);
varargout = cell(nout,1);
for i=1:nout
    if i > length(x)
        varargout(i) = {1};
    else
        varargout(i) = {x(i)}; 
    end
end

end

