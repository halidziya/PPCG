% clear;
load ..\data\crism\crism_mineral_data_wbackground.mat
Xorg = X;
%X=igmm_normalize(X(Y==4,:),50);
%Y=Y(Y==4);
%X=igmm_normalize(XX(YY==1,:),20);
%Y=YY(YY==1);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'crism');
igmm_mkdir(folder);
prefix = strcat(folder,'/','crism');


num_sweeps = '200';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ',prefix ,''];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=0.01;
gam=1;
s=1;
Psi=(m-d-1)*eye(d)/s;

igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
cla
scatter(X(:,1),X(:,2),40,labels+1,'.')

dpgmmres=evaluationTable(Y(Y~=0),labels(Y~=0))

for j=1:(max(labels)+1)
    if (table(j).npoints > 3)
    sigma = table(j).cholsigma'*table(j).cholsigma;
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],1.9,0.5);
    plot_gaussian_ellipsoid(table(j).mu(1:2),table(j).scatter(1:2,1:2)/table(j).npoints,'-',[0 0 0],1.8,2);
    end
end









% %% Experiments
% 
% hold off;
% for i=1:length(table)
% X2(:,i)=mvnpdf(X,table(i).mu',table(i).cholsigma*table(i).cholsigma')*table(i).npoints;
% end
% at=squareform(pdist(X2','cosine'))
% at=(cov(X2));
% at=at-min(min(at));
% at=at+0.01;
% %at=1./at;
% nodes=kamada_kawai_spring_layout(sparse(at))
% scatter(nodes(:,1),nodes(:,2));
% for i=1:size(X2,1)
% X3(i,:)=(X2(i,:)/sum(X2(i,:)))*nodes+mvnrnd(zeros(1,2),0.0000001*eye(2));
% end
% scatter(X3(:,1),X3(:,2),40,labels)