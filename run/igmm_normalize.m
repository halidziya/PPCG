function X=igmm_normalize(X,ndims)
    if (~exist('ndims','var'))        
        ndims = size(X,2);
    end
    if (size(X,2)<ndims)
        ndims = size(X,2);
    end
    
    X=X-ones(size(X,1),1)*mean(X,1);
    cc=cov(X);
    [vv dd]=eig(cc);
    X=X*vv(:,end:-1:(end-(ndims-1))); %d=20
    sg=std(X,[],1);
    X=X./(ones(size(X,1),1)*sg);
end