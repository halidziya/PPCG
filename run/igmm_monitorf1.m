%%Exeriment Folder
experimentname = 'toy';
prefix = ['experiments/' experimentname '/'];
if ~(exist(prefix, 'dir') == 7)
    mkdir(prefix);
end

addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian
addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian\include
addpath C:\Users\hzyereba\Desktop\JChang\JChang\common


macf1=[];
numtables=[];
likelihoods=[];
%% Generation
for iter=1:1
for D=2:4;
N1 = 1000;
N2 = 10000;
N3 = 100;
S  = 1;
d1 = mvnrnd(zeros(1,D),eye(D,D)/S,N1);
d2 = mvnrnd(10*ones(1,D)/sqrt(D),eye(D,D)/S,N2);
d3 = mvnrnd(-15*ones(1,D)/sqrt(D),eye(D,D)/S,N3);
X = [d1;d2;d3];
Y = [ones(N1,1);2*ones(N2,1);3*ones(N3,1)];
subplot(2,2,1);
scatter(X(:,1),X(:,2),4,Y);


%% Inference 
Psi = cov(X,1);
mu0 = mean(X,1);
m   = D+3;
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
likelihoods(D,method,iter,:) = readMyMat([prefix '_igmm.likelihood']);
subplot(2,2,method+1);
scatter(X(:,1),X(:,2),40,labels,'.')
title([ 'PPCG: ' num2str(macf1(D,method,iter))]);



%Call
slice_createBinaryFiles([prefix '/toy'],X,Psi,mu0,m,k0,gamma);
cmd = ['dpsl.exe ',data,' ',meanp,' ',psip,' ',params,' ',NITER,' ',BURNIN,' ',NSAMPLE];
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
subplot(2,2,method+1);
scatter(X(:,1),X(:,2),40,labels,'.')
title([ 'DPSL: ' num2str(macf1(D,method,iter))]);





[labels,E]=run_dpgmm_subclusters(X', 10, false, 8, false, false, 1, 500, 500);
method = 3;
labels = align_labels(labels);
f1s=evaluationTable(Y,labels);
macf1(D,method,iter)=table2array(f1s(1,1));
numtables(D,method,iter) = length(unique(labels));
likelihoods(D,method,iter,:)=E;
subplot(2,2,method+1);
scatter(X(:,1),X(:,2),40,labels,'.')
title([ 'JChang Sampler: ' num2str(macf1(D,method,iter))]);
   
drawnow;
end
end