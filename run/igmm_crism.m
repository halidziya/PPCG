clear;
load ..\data\crism\crism_mineral_data.mat
Xorg = X;
X=igmm_normalize(X,20);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'crism');
igmm_mkdir(folder);
prefix = strcat(folder,'/','crism');


num_sweeps = '200';
maxtables  =  '50';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ', maxtables , ' ',prefix];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=1;
gam=1;
s=1;
Psi=(m-d-1)*eye(d)/s;
igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
scatter(X(:,1),X(:,2),40,labels,'.')

dpgmmres=evaluationTable(Y(Y~=0),labels(Y~=0))

for j=1:(max(labels)+1)
    if (table(j).npoints > 25)
    sigma = table(j).cholsigma*table(j).cholsigma';
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],2,1);
    end
end