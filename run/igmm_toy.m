clear;
load ..\data\toy\ToyData_I2GMM_journal_3rare_3normal_3D_final.mat
%Xorg = X;
%X=[1.01 1.01;0.99 0.99; 1 1; 0 0 ;  -1.01 -1.01; -1 -1;-0.99 -0.99;1 -1];
%Y = [1;1;1;2;3;3;3;4];
X=igmm_normalize(X);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'toy');
igmm_mkdir(folder);
prefix = strcat(folder,'/','toy');


num_sweeps = '2500';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps  , ' ',prefix];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=0.01;
gam=1;
s=d^(1/2);
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
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',cc(j,:)/2,1,2);
    end
    
end
%axis([-2 2 -2 2])