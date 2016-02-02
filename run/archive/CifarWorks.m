load data\cifar-10-batches-mat\test_batch.mat;
X=double(data);
X=X-ones(size(X,1),1)*mean(X,1);
cc=cov(X);
[vv dd]=eig(cc);
novar=sum(diag(dd)<eps);
X=X*vv(:,end:-1:(novar+1));
sg=std(X,[],1);
X=X./(ones(size(X,1),1)*sg);

patch = reshape((1:3072),32,32,3);
patch = patch(10:17,10:17,1);
patch = patch(:);
X=X(randi(10000,5000,1),patch);
X=2*((X-repmat(min(X),size(X,1),1))./repmat(max(X)-min(X),size(X,1),1)-0.5);
lbs=kmeans(X,20);

options = [];
options.Fisherface = 1;
[eigvector, eigvalue] = LDA(lbs, options, X);
X2 = X*eigvector;


for i=1:4
    for j=5:8
        subplot(4,4,(i-1)*4+(j-4));
        %plot(X(:,at(i)),X(:,at(j)),'.')
         scatter(X2(:,i),X2(:,j),50,lbs,'.');
    end
end