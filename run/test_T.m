k0=0.01;
gam=10;
s=2;
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
    

scatter(X(:,1),X(:,2),10,labels+2,'.')

dots = (meshgrid(-2.5:0.05:2.5,-2.5:0.05:2.5));
dits = dots';
dots = [dits(:) dots(:)];

likelihoods = zeros(length(tables)+1,size(dots,1));
for j=1:(length(tables))
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
