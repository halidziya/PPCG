clear;
load data\toy\ToyData_I2GMM_journal_3rare_3normal_3D_final.mat
Xorg = X;
X=igmm_normalize(X);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'toy');
igmm_mkdir(folder);
prefix = strcat(folder,'/','toy');


num_sweeps = '1000';
burn_in = '100'; %Not active now
step = '100'; %Not active now
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps,' ', burn_in,' ',prefix,' ',step];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=0.5;
gam=0.01;
s=1;
Psi=(m-d-1)*eye(d)/s;
igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
scatter(Xorg(:,1),Xorg(:,2),40,labels,'.')