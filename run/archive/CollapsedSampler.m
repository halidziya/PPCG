function labels=CollapsedSampler(X,s,alpha,kappa) %I didn't add other params yet , mu0 is zero
    n = size(X,1);
    d = size(X,2);
    Psi = eye(d)/s;
    labels = kmeans(X,max(1,floor(size(X,1)/5))); % Initial 5 points in each class
    lik0 = mymvtpdf(X,s);
    lik0
    if (n>1)
        
        %Initial
        ntable = length(unique(labels));
        ns = zeros(1,ntable);
        xbar = zeros(ntable,d);
        ss   = zeros(d,d,ntable);
        
        for i=1:ntable
            idx = find(labels==i);
            ns(i) = length(idx);    
            xbar(i,:) = mean(X(idx,:),1);
            if (length(idx)>1)
                ss(:,:,i) = n*cov(X(idx,:)) + (ns(i)*kappa)*xbar(i,:)'*xbar(i,:)/(ns(i)+kappa);
            else
                ss(:,:,i) = zeros(d,d);
            end
        end
        
        
        % Sample
        for iter=1:5
            for i=1:n
                
                % Remove point
                y=labels(i);
                diff = X(i,:) - xbar(y,:);
                xbar(labels(i),:) = (xbar(y,:)*ns(y) - X(i,:))/(ns(y)-1);
                ss(:,:,y) = ss(:,:,y) - diff'*diff/(ns(y)/(ns(y)-1));
                ns(y)=ns(y)-1;
                
                
                
                % Calculate Likelihood
                liketable = zeros(1,ntable+1);
                for j=1:ntable
                    if (ns(j)==0)
                        liketable(j)=0;
                    else
                    liketable(j)=mvnpdf(X(i,:),xbar(j,:)*(ns(j)/(ns(j)+kappa)),(Psi + ss(:,:,j))*((ns(j)+kappa+1)/((kappa+ns(j))*(kappa + n + d - 1))) )*ns(j);
                    end
                end
                liketable(ntable+1)=lik0(i)*alpha;
                
                %Sample
                liketable
                liketable = liketable / sum(liketable);
                y = find(mnrnd(1,liketable)); %Selected
                try
                labels(i) = y;
                catch
                    keyboard
                end
                %Add point
                if (y > ntable)
                    ntable = ntable+1;
                    ns(ntable)=0;
                    xbar(ntable,:) = zeros(1,d);
                    ss(:,:,ntable) = zeros(d,d);
                end
                
                ns(y)=ns(y)+1;
                xbar(labels(i),:) = (xbar(y,:)*(ns(y)-1) + X(i,:))/(ns(y)); 
                diff = X(i,:) - xbar(y,:);
                ss(:,:,y) = ss(:,:,y) + diff'*diff/((ns(y))/(ns(y)-1));
                             
            end
        end
        
    NK =length(unique(labels));
    used=1:max(labels);
    used(unique(labels))=1:NK;
    labels=used(labels);
    end
end