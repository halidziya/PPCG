function likelihoods=igmm_hierarchylikelihood(X,model)
    k0=0.01;
    d= size(X,2);
    n= size(X,1);
    mu0 = mean(X);
    s=1;
    m = d+3;
    Psi=eye(d)/s;
    likelihoods=zeros(length(model),size(X,1));
    for i=1:length(model)
        mus = model{i}{1};
        covs = model{i}{2};
        np   = model{i}{3};
        postmus =  (mus .* repmat(np,1,d) + k0*repmat(mu0,length(np),1)) ./ (repmat(np,1,d)  + k0);
        mixture = zeros(length(np),n);
        for j=1:size(mus,1)
        postcov = (Psi + squeeze(covs(j,:,:))*(np(j)-1) + (mus(j,:)-mu0)'*(mus(j,:)-mu0)*((k0*np(j))/(k0+np(j))))/(np(j)+1);
        mixture(j,:) = mvnpdf(X,postmus(j,:),postcov)*np(j);
        display(j)
        end
        likelihoods(i,:) = sum(mixture,1); 
    end
end