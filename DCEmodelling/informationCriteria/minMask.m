function [ mask ] = minMask( varargin )
%function [ mask ] = minMask( varargin )
%   
%   elementwise arg min for masks
%   usage:
%       mask = minMask( model1, model2, model3 );
%   mask will contain the index of the model with smallest value for each
%   voxel
%
%   Yannick Bliesener 2018

dims = size( varargin{1} );

for k=1:nargin
    varargin{k} = reshape(varargin{k}, [], 1);
end

mask    = ones(size(varargin{1}));
curMin  = varargin{1};

for k=2:nargin
    mask(curMin > varargin{k}) = k;
    curMin = min(curMin, varargin{2});   
end

mask = reshape(mask, dims);

end

