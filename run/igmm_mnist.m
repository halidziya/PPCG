load('C:\Users\halidziya\Desktop\I2GMM\UncollapsedSampler\run\data\mnist\mnist_all.mat')
X=[train0; train1; train2; train3; train4; train5; train6; train7; train8; train9];
Y=[ones(size(train0,1),1); 2*ones(size(train1,1),1); 3*ones(size(train2,1),1); 4*ones(size(train3,1),1); 5*ones(size(train4,1),1); 6*ones(size(train5,1),1); 7*ones(size(train6,1),1); 8*ones(size(train7,1),1); 9*ones(size(train8,1),1); 10*ones(size(train9,1),1) ];

Xorg = X;
X=igmm_normalize(double(X),10);
%Y=Y(Y==4);
%X=igmm_normalize(XX(YY==1,:),20);
%Y=YY(YY==1);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'mnist');
igmm_mkdir(folder);
prefix = strcat(folder,'/','mnist');


num_sweeps = '200';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ',prefix ,''];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+2;
mu0 = mean(X);
k0=0.1;
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

%5 dim elapsed time 42s F1 0.51 0.52
%50 dim elapsed time 20 min F1 0.18
%50 dim elapsed time 249 min F1 0.18 ==> Split merge
%5 dim elapsed time 76s F1 0.52 ==> Split merge
%10 dim elapsed time 87s F1 0.61 - 0.56 ==> heterosampler
%10 dim elapsed time 17min F1 0.49421 ==> Split merge

