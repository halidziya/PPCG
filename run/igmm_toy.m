clear;
load ..\data\toy\ToyData_I2GMM_journal_3rare_3normal_3D_final.mat
Xorg = X;
X=igmm_normalize(X);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'toy');
igmm_mkdir(folder);
prefix = strcat(folder,'/','toy');


num_sweeps = '100';
maxtables  =  '10';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ', maxtables , ' ',prefix];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = 2*d+2;
mu0 = mean(X);
k0=0.5;
gam=0.1;
s=10;
Psi=(m-d-1)*eye(d)/s;
igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
clf
scatter(X(:,1),X(:,2),40,labels,'.')


for j=1:(max(labels)+1)
    if (table(j).npoints > 0)
    sigma = table(j).cholsigma*table(j).cholsigma';
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],2,1);
    end
end