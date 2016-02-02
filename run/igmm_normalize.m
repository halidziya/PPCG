function X=igmm_normalize(X)
    X=X-ones(size(X,1),1)*mean(X,1);
    cc=cov(X);
    [vv dd]=eig(cc);
    X=X*vv;
    sg=std(X,[],1);
    X=X./(ones(size(X,1),1)*sg);
end