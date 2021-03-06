experiments='experiments/';
folder = strcat(experiments,'all2');
igmm_mkdir(folder);
[files names] =  igmm_datasets('..\data'); % Traverse in folder
MAXITER=1;
elapsed_time = zeros(length(files),5,MAXITER);
macf1        = zeros(length(files),5,MAXITER);
micf1        = zeros(length(files),5,MAXITER);

addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian
addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian\include
addpath C:\Users\hzyereba\Desktop\JChang\JChang\common


for datai=1:length(names)

    prefix = char(strcat(folder,'/',names(datai)));
    mkdir([prefix,'\plots\']);
    run(files{datai});
    X=igmm_normalize(X,20);
    

    num_sweeps = '2000';
    data=[prefix,'.matrix'];
    prior=[prefix,'_prior.matrix'];
    params=[prefix,'_params.matrix'];
    
    fprintf(1,'\nIGMM is running...\n');
    d=size(X,2);
    m = d+3;
    mu0 = mean(X);
    k0=0.01;
    gam=1;
    s=1;
    Psi=eye(d)/s; %(m-d-1);
    for iter=1:MAXITER
    igmm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam);
    cmd = ['igmm.exe ',data,' ',prior,' ',params,' ',num_sweeps  , ' ',prefix];
    tic;
    system(cmd);
    elapsed_time(datai,1,iter)=toc;
    
    [table labels]=igmm_readOutput([prefix '_igmm.rest']);
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    
    macf1(datai,1,iter)=table2array(f1s(1,1));
    micf1(datai,1,iter)=table2array(f1s(1,2));
   
    
    cmd = ['igmm_kd.exe ',data,' ',prior,' ',params,' ',num_sweeps  , ' ',prefix];
    tic;
    system(cmd);
    elapsed_time(datai,2,iter)=toc;
    
    [table labels]=igmm_readOutput([prefix '_igmm.rest']);
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    
    macf1(datai,2,iter)=table2array(f1s(1,1));
    micf1(datai,2,iter)=table2array(f1s(1,2));
    
    subplot(2,2,1);
    cla
    scatter(X(:,1),X(:,2),40,labels+2,'.')
    for j=1:(max(labels)+1)
        if (table(j).npoints > 2)
            sigma = table(j).cholsigma'*table(j).cholsigma;
            plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],2,0.5);
            plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],1,1);
        end
    end
    title([ 'IGMM HeteroCollapsed Sampler: ' num2str(macf1(datai,1,iter))]);
    
    burn_in = '500';
    step = '100';
    dpm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam,1);
    cmd = ['dpm64.exe ',data,' ',prior,' ',params,' ',num_sweeps,' ', burn_in,' ',prefix,' ',step];
    fprintf(1,'\nDPGMM is running...\n');
    tic;
    system(cmd);
    elapsed_time(datai,3,iter)=toc;

    [tables customers klabels]=dpm_readOutput(prefix);
    labels = klabels(:,end);
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    macf1(datai,3,iter)=table2array(f1s(1,1));
    micf1(datai,3,iter)=table2array(f1s(1,2));
    
    subplot(2,2,2);
    scatter(X(:,1),X(:,2),40,labels,'.')
    title([ 'Collapsed Sampler: ' num2str(macf1(datai,3,iter))]);
    
    
    
    
    
    tic;
    labels=run_dpgmm_subclusters(X', 1, false, 8, false, false, 1, 2000, 2000);
    elapsed_time(datai,4,iter)=toc;
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    macf1(datai,4,iter)=table2array(f1s(1,1));
    micf1(datai,4,iter)=table2array(f1s(1,2));
    subplot(2,2,3);
    scatter(X(:,1),X(:,2),40,labels,'.')
    title([ 'JChang Sampler: ' num2str(macf1(datai,4,iter))]);
    
    
       
    tic;
    hierlabels = igmm_hierarchy(X);
    elapsed_time(datai,5,iter)=toc;
    labels = igmm_hierarchylabels(X,hierlabels);
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    macf1(datai,5,iter)=table2array(f1s(1,1));
    micf1(datai,5,iter)=table2array(f1s(1,2));
    subplot(2,2,4);
    scatter(X(:,1),X(:,2),40,labels,'.')
    title([ 'Hier Sampler: ' num2str(macf1(datai,5,iter))]);
    
    
    
    print(cell2mat(strcat(prefix,'\plots\',names(datai),'_',num2str(iter))),'-dpng','-r300');
    end
end