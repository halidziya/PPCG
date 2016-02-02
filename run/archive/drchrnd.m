%Sample from Dirichlet Distribution
%r = drchrnd(eta,dim)
function r = drchrnd(eta,dim)
    if size(eta,2) == 1 % Symmetric dirichlet
        r=gamrnd(eta*ones(1,dim),1);
    else
        r=gamrnd(eta,1);
    end
    % Always sum to 1
    r = r/sum(r);
    % Abs value is used to avoid -0 problem
    r(end) = abs(1 - sum(r(1:end-1)));
end
