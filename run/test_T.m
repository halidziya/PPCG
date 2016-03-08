k0=0.1;
gam=100;
s=10;
Psi=(m-d-1)*eye(d)/s;
m=d+2;

cla

for j=1:(max(labels)+1)
    if (table(j).npoints > 2)
        sigma = table(j).cholsigma'*table(j).cholsigma;
        plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],2,0.5);
        plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],1,1);
    end
end
    

scatter(X(:,1),X(:,2),40,labels+2,'.')

dots = (meshgrid(-2:0.02:2,-2:0.02:2));
dits = dots';
dots = [dits(:) dots(:)];

likelihoods = zeros(length(table)+1,size(dots,1));
for j=1:(length(table))
        if (table(j).npoints > 2)
            sigma = table(j).cholsigma'*table(j).cholsigma;
            likelihoods(j,:) = mvnpdf(dots,table(j).mu(1:2)',sigma(1:2,1:2))*table(j).npoints;
        end
end
tC = Psi*((k0 + 1) / ((k0)*(m-d+1)));
likelihoods(j+1,:) = mvtpdf(dots,tC(1:2,1:2),m - d + 1)*gam;
[a b]=max(likelihoods)

scatter(dots(:,1),dots(:,2),150,b,'.');
igmm_colorSettings
