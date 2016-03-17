experiments='experiments/';
folder = strcat(experiments,'parallel');
igmm_mkdir(folder);
[files names] =  igmm_datasets('..\data'); % Traverse in folder
MAXITER=2;
elapsed_time = zeros(length(files),2,MAXITER);
macf1        = zeros(length(files),2,MAXITER);
micf1        = zeros(length(files),2,MAXITER);
numtables        = zeros(length(files),2,MAXITER);

addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian
addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian\include
addpath C:\Users\hzyereba\Desktop\JChang\JChang\common


for datai=1:length(names)

    prefix = char(strcat(folder,'/',names(datai)));
    mkdir([prefix,'\plots\']);
    run(files{datai});
    X=igmm_normalize(X,50);
    
    num_sweeps = '2000';
    data=[prefix,'.matrix'];
    prior=[prefix,'_prior.matrix'];
    params=[prefix,'_params.matrix'];
    
    fprintf(1,'\nIGMM is running...\n');
    d=size(X,2);
    m = d+3;
    mu0 = mean(X);
    k0=1;
    gam=1;
    %s=1/d;
    Psi=eye(d)*m;
    for iter=1:MAXITER
    igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);
    cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps  , ' ',prefix];
    tic;
    system(cmd);
    elapsed_time(datai,1,iter)=toc;
    fprintf(1,'Reading...\n');
    [table,labels]=igmm_readOutput([prefix '_igmm.rest']);
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    
    macf1(datai,1,iter)=table2array(f1s(1,1));
    micf1(datai,1,iter)=table2array(f1s(1,2));
   numtables(datai,1,iter) = length(table);
   
    subplot(1,2,1);
    cla;
    scatter(X(:,1),X(:,2),40,labels,'.')
        for j=1:(max(labels)+1)
        if (table(j).npoints > 2)
            sigma = table(j).cholsigma'*table(j).cholsigma;
            plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],2,0.5);
            plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],1,1);
        end
        end
    
    title([ 'IGMM HeteroCollapsed Sampler: ' num2str(macf1(datai,1,iter))]);
    
    
    tic;
    labels=run_dpgmm_subclusters(X', 1, false, 8, false, false, 1, 2000, 2000);
    %labels=run_dpgmm_fsd(X',1,false,8,1,40,40);
    elapsed_time(datai,2,iter)=toc;
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    macf1(datai,2,iter)=table2array(f1s(1,1));
    micf1(datai,2,iter)=table2array(f1s(1,2));
    numtables(datai,2,iter) = length(unique(labels));
    subplot(1,2,2);
    scatter(X(:,1),X(:,2),40,labels,'.')
    title([ 'JChang Sampler: ' num2str(macf1(datai,2,iter))]);
    
       
    print(cell2mat(strcat(prefix,'\plots\',names(datai),'_',num2str(iter))),'-dpng','-r300');
    end
end