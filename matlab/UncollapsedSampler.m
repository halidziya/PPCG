function [labels M]=UncollapsedSampler(X,Y)
%load('ToyData_I2GMM_journal_3rare_3normal_3D_final.mat')
%parpool(7);

%Whiten the Data
X=X-ones(size(X,1),1)*mean(X,1);
cc=cov(X);
[vv dd]=eig(cc);
novar=sum(diag(dd)<eps);
X=X*vv(:,end:-1:(novar+1));
sg=std(X,[],1);
X=X./(ones(size(X,1),1)*sg);
%X=X(:,1:2); %2D for now


%Plot
subplot(2,1,1);
scatter(X(:,1),X(:,2),50,Y,'.');
title('Original')
%X=X(Y<5&Y>1,:);

% Initial KMEANS
NP = size(X,1);
NK = 4;
labels = kmeans(X,NK);

% Parameters
d=size(X,2);
m = d+2;
s=(100)^(1/d);
kappa0=1/s;
alpha = 1;

% Preallocation of Matrices
mus = zeros(NK,d);
sigmas =  zeros(d,d,NK);

prevNK=NK;
% New class does not change prob , due to fixed prior
newclass = mymvtpdf(X,s);  
for iter=1:100
    fprintf(1,'Iter: %d NK: %d\n',iter,NK);
    
    % Update statistics according to posterior
    ns = zeros(NK,1);
     for i=1:NK
        idx = find(labels==i);
        n = length(idx);
        ns(i) = n;
        xbar = mean(X(idx,:),1);
        if (n>1)
            sigmas(:,:,i) = iwishrnd(n*cov(X(idx,:))+eye(d)/s+((kappa0*n)/(kappa0+n))*(xbar'*xbar),m+n);
        else
            fprintf('Here %d\n',n);
            sigmas(:,:,i) = iwishrnd(eye(d)/s+((kappa0*n)/(kappa0+n))*(xbar'*xbar),m+n);
        end
        
        
        if (n>0)
            mus(i,:) = mvnrnd(xbar*n/(n+kappa0),sigmas(:,:,i)/(kappa0+n));
        else
            mus(i,:) = mvnrnd(zeros(d,1),sigmas(:,:,i)/(kappa0));
        end
        
     end
    subplot(1,1,1);
    scatter(X(:,1),X(:,2),50,labels,'.');
    
    hold on;
    scatter(mus(1:NK,1),mus(1:NK,2),100,'rx')
    %if (sum(labels(labels>prevNK))>0)
    %scatter(X(labels>prevNK,1),X(labels>prevNK,2),850,labels(labels>prevNK),'.');
    %end
    plot_gaussian_ellipsoid([0,0],eye(2),':',0.1*ones(3,1),2,0.5)
    for i=1:NK
        if (ns(i) > 100)
        plot_gaussian_ellipsoid(mus(i,1:2),sigmas(1:2,1:2,i),':',0.5*ones(3,1),2,2)
        end
    plot_gaussian_ellipsoid(mus(i,1:2),sigmas(1:2,1:2,i),':',0.2*ones(3,1),log(ns(i))/2,2)
    end
    hold off;
    axis([-4 4 -4 4])
    drawnow
    

      
      
      
     
    % Calculate Likelihoods
    liketable=zeros(NP,NK);
    counts = ns;
    counts'
    ns=drchrnd([ns ;alpha],NK+1);
    for j=1:NK
    %Slicing for MATLAB
    
    mu=mus(j,:);
    sigma=sigmas(:,:,j);
    n = ns(j);
    
    % Data-Parallel Calculation
    parfor i=1:NP
        liketable(i,j)=mvnpdf(X(i,:),mu,sigma)*(n);
    end
    end

    % Sample
    probs=[ liketable newclass.*ns(end) ];
    probs=probs./repmat(sum(probs,2),1,size(probs,2));
    [dum,idx]=sort(mnrnd(1,probs),2,'descend');
    labels = idx(:,1);
    if (sum(labels==(NK+1)) > 0)
        labels(labels==(NK+1)) = NK+CollapsedSampler(X(labels==(NK+1),:),s,alpha,kappa0);
    end
    prevNK=NK;
%     % New Tables
%     ntables = find(labels==NK+1);
%     labels(ntables) = NK+(1:length(ntables));
%     NK = NK+length(ntables);
%     
    %Remove Tables

    NK =length(unique(labels));
    used=1:max(labels);
    used(unique(labels))=1:NK;
    labels=used(labels);

      drawnow
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if iter == 1;
          imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf);
      else 
          if(mod(iter,5)==0)
          imwrite(imind,cm,'test.gif','gif','WriteMode','append');
          end
      end
    
    M(iter) = getframe;
end
end



