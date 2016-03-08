function [coords ww bias]=visSVM(X,Y,selected)
at=fitcsvm(X,Y==selected);
w=sum(at.SupportVectors.*repmat(at.SupportVectorLabels.*at.Alpha,1,size(X,2)));
bias = at.Bias;
old=w;
w=w/norm(w);
Xremain = (X-X*(w'*w));
[a ~]=eigs(cov(Xremain));
ww=[old' a(:,1)];
coords = X*ww;
coords(:,1)=(coords(:,1)+bias);
cla;
hold on;
scatter(coords(:,1),coords(:,2),40,Y,'.');
scatter(coords(Y==selected,1),coords(Y==selected,2),40,Y(Y==selected),'o');
hold off;
end