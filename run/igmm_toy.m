clear;
%load ..\data\toy\Toy.mat
%Xorg = X;
%X=[-1 -1; -0.99 -0.99; -1.01 -1.01; -1 -0.99;1 1; 0.99 0.99; 1.01 1.01; 1 0.99;-2 2;-2.01 2;];
%X=X;
%Y = [1 1 1 1 2 2 2 2 3];
%load fisheriris;
%X=meas;
%X=igmm_normalize(X);
d=2;
X=[mvnrnd(ones(1,d),eye(d)/20,1000);mvnrnd(-0.1*ones(1,d),eye(d)/20,100);mvnrnd(-1.1*ones(1,d),eye(d)/20,50)];
%X=igmm_normalize(X);
experiments='experiments/';
folder = strcat(experiments,'toy');
igmm_mkdir(folder);
prefix = strcat(folder,'/','toy');


num_sweeps = '100.';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps  , ' ',prefix];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+3;
mu0 = mean(X);
k0=1;
gam=1;
s=1;
Psi=eye(d)/s;
igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);

tic;
system(cmd);
elapsed_time(1)=toc;

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
clf
scatter(X(:,1),X(:,2),40,floor(labels)+1,'.')

igmm_colorSettings;
for j=1:(max(labels)+1)
    if (table(j).npoints > 1)
    sigma = table(j).cholsigma'*table(j).cholsigma;
    %plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],log(table(j).npoints)/5,2);
    plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',cc(j,:)/2,1,2);
    at = (cov(X(labels==(j-1),1:2))+0.001*eye(2));
    bet = mean(X(labels==(j-1),1:2));
    plot_gaussian_ellipsoid(bet,at,'-',cc(j,:),1.1,1.5);
    for i=1:10
        nn = sum(labels==(j-1));
        sampleCov=iwishrnd(at*nn+Psi(1:2,1:2)*m,nn+m);
        sampleMean=mvnrnd((bet*nn + mu0(1:2)*k0)/(k0+nn),sampleCov/(k0+nn));
        plot_gaussian_ellipsoid(sampleMean,sampleCov,'-',[0.7 0.7 0.7],0.5,0.5);
        plot_gaussian_ellipsoid(sampleMean,sampleCov,'-',[0.7 0.7 0.7],0.5,0.5);
    end
    end
end

%axis([-2 2 -2 2])