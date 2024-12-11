function dec = sample_dec(T,varargin)
sz = zeros(nargin,1);
sz(1) = T;
for n = 2:nargin
    sz(n) = varargin{n-1};
end
dec = rand(sz');
end