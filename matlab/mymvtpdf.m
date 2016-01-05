function y = mymvtpdf(X,s)
%A slight modification of MATLAB student-t 
if (~exist('s'))
    s=10;
end
d=size(X,2);
df = 3;
C=eye(d)*((s+1)/s);


if nargin<1
    error(message('mymvtpdf:TooFewInputs'));
elseif ndims(X)~=2
    error(message('mymvtpdf:InvalidData'));
end

% Get size of data.  Column vectors provisionally interpreted as multiple scalar data.
[n,d] = size(X);
if d<1
    error(message('stats:mvtpdf:TooFewDimensions'));
end

% Special case: try to interpret X as a row vector if it was a column.
if isvector(X) && (size(C,1) == n)
    X = X';
    [n,d] = size(X);
end

sz = size(C);
if sz(1) ~= sz(2)
    error(message('stats:mvtpdf:BadCorrelationNotSquare'));
elseif ~isequal(sz, [d d])
    error(message('stats:mvtpdf:InputSizeMismatchC'));
end

% Make sure C is a valid covariance matrix
[R,err] = cholcov(C,0);
if err ~= 0
    error(message('stats:mvtpdf:BadCorrelationSymPos'));
end

if ~(isscalar(df) || (isvector(df) && length(df) == n))
    error(message('stats:mvtpdf:InputSizeMismatchDF'));
elseif any(df <= 0)
    error(message('stats:mvtpdf:InvalidDF'));
end
df = df(:);

% Create array of standardized data, and compute log(sqrt(det(Sigma)))
Z = X / R;
logSqrtDetC = sum(log(diag(R)));

logNumer = -((df+d)/2) .* log(1+sum(Z.^2, 2)./df);
logDenom = logSqrtDetC + (d/2)*log(df*pi);
y = exp(gammaln((df+d)/2) - gammaln(df/2) + logNumer - logDenom);
