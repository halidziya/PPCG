%Xorg = X;
X=igmm_normalize(X,20);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'crism');
igmm_mkdir(folder);
prefix = strcat(folder,'/','crism');
elapsed_time=zeros(2,1)

num_sweeps = '2000';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ',prefix];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+20;
mu0 = mean(X);
k0=0.0001;
gam=10^4;
s=10^3;
Psi=(m-d-1)*eye(d)/s;
igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
cla
scatter(X(:,1),X(:,2),40,labels+2,'.')

dpgmmres=evaluationTable(Y(Y~=0),labels(Y~=0))

for j=1:(max(labels)+1)
    if (table(j).npoints > 1)
    sigma = table(j).cholsigma'*table(j).cholsigma;
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],2,0.5);
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],1,2);
    end
end