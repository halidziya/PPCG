[X Y]=meshgrid(-2:0.01:2,-2:0.01:2);X=[X(:) Y(:)];
tval=mymvtpdf(X,eye(2),1);
C1 = iwishrnd(eye(2)/10,4);
C2 = iwishrnd(eye(2)/10,4);
n1val=mvnpdf(X,[-1,0],C1);
n2val=mvnpdf(X,[1,0],C2);
fraction = tval./(200*n1val+2*n2val+tval);
cla;
contourf(-2:0.01:2,-2:0.01:2,reshape(fraction,401,401))
fraction(201*401+201)
X=mvnrnd([1,0],C2,100);
hold on;plot(X(:,1),X(:,2),'r.')
X=mvnrnd([-1,0],C1,100);
hold on;plot(X(:,1),X(:,2),'g.')