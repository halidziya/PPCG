function labels=igmm_uncollapsed(X)

experiments='experiments/';
folder = strcat(experiments,'fun');
igmm_mkdir(folder);
prefix = strcat(folder,'/','fun');


num_sweeps = '2000';
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

end