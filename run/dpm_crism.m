clear;
load ..\data\crism\crism_mineral_data.mat
Xorg = X;
X=igmm_normalize(X,20);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'crism');
igmm_mkdir(folder);
prefix = strcat(folder,'/','crism');

num_sweeps = '2000';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
fprintf(1,'\nIGMM_collapsed is running...\n');

burn_in = '1000';
d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=0.001;
gam=1;
s=1;
step = '100';
Psi=(m-d-1)*eye(d)/s;
dpm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,1,gam);
cmd = ['dpm64.exe ',data,' ',prior,' ',params,' ',num_sweeps,' ', burn_in,' ',prefix,' ',step];
fprintf(1,'\nDPGMM is running...\n');
tic;
system(cmd);
elapsed_time=toc;

[tables customers klabels]=dpm_readOutput(prefix);
labels = klabels(:,end);
scatter(X(:,1),X(:,2),40,labels,'.')

dpgmmres=evaluationTable(Y(Y~=0),labels(Y~=0))