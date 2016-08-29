function best=align_labels(labels)
if ~isempty(find(labels==0,1))
    labels = labels+1; % Matlab index
end

    n = size(labels,1); % Each column is a sampled labels
    nsample = size(labels,2);
    labs=[];
    for i=2:nsample % First One is the reference
        labs(unique(labels(:,i)))=1:length(unique(labels(:,i))); % If it is not ordered
        cfm=confusionmat(labels(:,1),labels(:,i));
        cost = n - cfm;
        allignment = munkres(cost);
        [sorted inverseal]=sort(allignment);
        labels(:,i)=inverseal(labs(labels(:,i)));
    end
    
    best = mode(labels,2);
    
end