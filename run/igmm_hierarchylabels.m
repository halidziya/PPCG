function labels=igmm_hierarchylabels(X,hierlabels)
d=size(X,2);
ids = unique(hierlabels,'rows');
newlabels =ones(size(hierlabels,1),1);
for i=1:size(hierlabels,1)
    [a,idx]=ismember(hierlabels(i,:),ids,'rows');
    newlabels(i) = idx;
end

mus = zeros(max(newlabels),d);
covs = zeros(max(newlabels),d,d);
np = zeros(max(newlabels),1);
for i=1:max(newlabels)
    mus(i,:) = mean(X(newlabels==i,:));
    covs(i,:,:) = cov(X(newlabels==i,:));
    np(i)=sum(newlabels==i);
    model{i}={mus(i,:),covs(i,:,:),np(i)};
end
[vals,labels] = max(igmm_hierarchylikelihood(X,model));
end