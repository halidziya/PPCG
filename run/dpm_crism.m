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
maxtables  =  '50';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['dpm64.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ', maxtables , ' ',prefix];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=1;
gam=1;
alpha=1;
s=1;
Psi=(m-d-1)*eye(d)/s;
dpm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,alpha,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[tables customers klabels]=dpm_readOutput(prefix);
labels = klabels(:,end);
scatter(X(:,1),X(:,2),40,labels,'.')

dpgmmres=evaluationTable(Y(Y~=0),labels(Y~=0))