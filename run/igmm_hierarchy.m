function labels=igmm_hierarchy(Xrem)
ndim=10;
if (size(Xrem,2)<ndim)
    ndim = size(Xrem,2);
end
X=Xrem(:,1:ndim);
%X=igmm_normalize(Xrem(:,1:ndim));
%X=igmm_normalize(XX(YY==1,:),20);
%Y=YY(YY==1);
igmm_colorSettings;

experiments='experiments/';
folder = strcat(experiments,'hier');
igmm_mkdir(folder);
prefix = strcat(folder,'/','hier');


num_sweeps = '2000';
data=[prefix,'.matrix'];
prior=[prefix,'_prior.matrix'];
params=[prefix,'_params.matrix'];
cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps , ' ',prefix ,''];
fprintf(1,'\nIGMM is running...\n');


d=size(X,2);
m = d+3;
mu0 = mean(X);
k0=0.01;
gam=1;
s=1;
Psi=(m-d-1)*eye(d)/s;

igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);
system(cmd);

[table labels]=igmm_readOutput([prefix '_igmm.rest']);
labels=labels+1;
%cla
%scatter(X(:,1),X(:,2),40,labels+1,'.')


np=histc(labels,unique(labels));
deeplabels=ones(size(labels,1),ceil(size(Xrem,2)/ndim)-1);
for i=1:length(np)
    if (np(i)>2 && size(Xrem,2)>ndim)
        dlabels = igmm_hierarchy(Xrem(labels==i,(ndim+1):end));
        deeplabels(labels==i,:)=dlabels;
    end
end
labels = [labels deeplabels];
end