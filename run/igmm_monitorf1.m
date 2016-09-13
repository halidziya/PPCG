%%Exeriment Folder
experimentname = 'toy';
prefix = ['experiments/' experimentname '/'];
if ~(exist(prefix, 'dir') == 7)
    mkdir(prefix);
end

addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian
addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian\include
addpath C:\Users\hzyereba\Desktop\JChang\JChang\common
addpath C:\Users\halidziya\Desktop\I2GMM\JChang\Gaussian
addpath C:\Users\halidziya\Desktop\I2GMM\JChang\Gaussian\include
addpath C:\Users\halidziya\Desktop\I2GMM\JChang\common

macf1=[];
numtables=[];
likelihoods=[];
elapsedtime = [];
%% Generation
for iter=1:10;
for D=2:20;
N1 = 10;
N2 = 10000;
N3 = 500;
S  = 1;
d1 = mvnrnd(zeros(1,D),eye(D,D)/S,N1);
d2 = mvnrnd(10*ones(1,D)/sqrt(D),eye(D,D)/S,N2);
d3 = mvnrnd(-15*ones(1,D)/sqrt(D),eye(D,D)/S,N3);
X = [d1;d2;d3];
Y = [ones(N1,1);2*ones(N2,1);3*ones(N3,1)];
subplot(2,2,1);
scatter(X(:,1),X(:,2),4,Y);


%% Inference 
mu0 = mean(X,1);
m   = D+3;
Psi = eye(D)*2;
k0  = 1;
gamma = 1;

%File Names
data=[prefix,'toy.matrix'];
prior=[prefix,'toy_prior.matrix'];
params=[prefix,'toy_params.matrix'];
psip=[prefix,'toy_psi.matrix'];
params=[prefix,'toy_params.matrix'];
NITER = '1000';
BURNIN = '800';
NSAMPLE = '10';


%Call
igmm_createBinaryFiles([prefix 'toy'],X,Psi,mu0,m,k0,gamma);
cmd = ['ppcg.exe ',data,' ',prior,' ',params,' ',NITER  , ' ',prefix,'toy'];
tic;
system(cmd);
fprintf(1,'Reading...\n');
[table,llabels]=igmm_readOutput([prefix 'toy_igmm.rest']);


method = 1;
labels = align_labels(readMyMat([prefix 'toy_igmm.labels']));
f1s=evaluationTable(Y,labels);
macf1(D,method,iter)=table2array(f1s(1,1));
numtables(D,method,iter) = length(unique(labels));
likelihoods(D,method,iter,:) = readMyMat([prefix 'toy_igmm.likelihood']);
elapsedtime(D,method,iter)=toc;
subplot(2,2,method+1);
scatter(X(:,1),X(:,2),40,labels,'.')
title([ 'PPCG: ' num2str(macf1(D,method,iter))]);



%Call
slice_createBinaryFiles([prefix '/toy'],X,Psi,mu0,m,k0,gamma);
cmd = ['dpsl.exe ',data,' ',meanp,' ',psip,' ',params,' ',NITER,' ',BURNIN,' ',NSAMPLE];
tic;
fprintf(1,[cmd , '\n']);
system(cmd);

%Combine multiple label samples
prediction=readMat([data '.labels']);
method = 2;
labels = align_labels(prediction');
f1s=evaluationTable(Y,labels);
macf1(D,method,iter)=table2array(f1s(1,1));
numtables(D,method,iter) = length(unique(labels));
likelihoods(D,method,iter,:)=readMat([data '.likelihood']);
elapsedtime(D,method,iter)=toc;
subplot(2,2,method+1);
scatter(X(:,1),X(:,2),40,labels,'.')
title([ 'DPSL: ' num2str(macf1(D,method,iter))]);




tic;
[labels,E]=run_dpgmm_subclusters(X', 1, false, 8, false, false, 1, 1000, 1000);
method = 3;
labels = align_labels(labels);
f1s=evaluationTable(Y,labels);
macf1(D,method,iter)=table2array(f1s(1,1));
numtables(D,method,iter) = length(unique(labels));
likelihoods(D,method,iter,:)=E(2:end);
elapsedtime(D,method,iter)=toc;
subplot(2,2,method+1);
scatter(X(:,1),X(:,2),40,labels,'.')
title([ 'JChang Sampler: ' num2str(macf1(D,method,iter))]);
   
drawnow;
end
end

plot(squeeze(mean((likelihoods(2,:,:,:)),3))')