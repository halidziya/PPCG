experiments='experiments/';
folder = strcat(experiments,'parallel');
igmm_mkdir(folder);
[files names] =  igmm_datasets('..\..\data'); % Traverse in folder
MAXITER=10;
elapsed_time = zeros(length(files),4,MAXITER);
macf1        = zeros(length(files),4,MAXITER);
micf1        = zeros(length(files),4,MAXITER);
numtables        = zeros(length(files),4,MAXITER);
effectiven        = zeros(length(files),4,MAXITER);

addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian
addpath C:\Users\hzyereba\Desktop\JChang\JChang\Gaussian\include
addpath C:\Users\hzyereba\Desktop\JChang\JChang\common

colormap(distinguishable_colors(25));
likelihood=[];
slikelihood=[];
changlikelihood=[];
for datai=1:length(names)

    prefix = char(strcat(folder,'/',names(datai)));
    mkdir([prefix,'\plots\']);
    run(files{datai});
    X=igmm_normalize(X,50);
    
    num_sweeps = '500';
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
    cmd = ['ppcg.exe ',data,' ',prior,' ',params,' ',num_sweeps  , ' ',prefix];
    tic;
    system(cmd);
    elapsed_time(datai,1,iter)=toc;
    fprintf(1,'Reading...\n');
    [table,llabels]=igmm_readOutput([prefix '_igmm.rest']);
    labels = align_labels(readMyMat([prefix '_igmm.labels']));
    p=histc(labels,unique(labels))/length(labels);
    effectiven(datai,1,iter) = exp(-sum(p.*log(p)));
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    
    macf1(datai,1,iter)=table2array(f1s(1,1));
    micf1(datai,1,iter)=table2array(f1s(1,2));
    numtables(datai,1,iter) = length(table);
    likelihood(:,iter) = readMyMat([prefix '_igmm.likelihood']);
    
    subplot(2,2,1);
    cla;
    scatter(X(:,1),X(:,2),40,llabels,'.')
    
    title([ 'IGMM HeteroCollapsed Sampler: ' num2str(macf1(datai,1,iter))]);
    
    
    data=[prefix,'.matrix'];
    prior=[prefix,'_prior.matrix'];
    params=[prefix,'_params.matrix'];
    psip=[prefix,'_psi.matrix'];
    meanp=[prefix,'_mean.matrix'];
    params=[prefix,'_params.matrix'];
    NITER = '500';
    BURNIN = '300';
    NSAMPLE = '10';
    slice_createBinaryFiles(prefix,X,Psi,mu0,m,k0,1);
    cmd = ['dpsl.exe ',data,' ',meanp,' ',psip,' ',params,' ',NITER,' ',BURNIN,' ',NSAMPLE];
    tic;
    system(cmd);
    elapsed_time(datai,3,iter)=toc;
    fprintf(1,'Reading...\n');
    labels = align_labels(readMyMat([data '.labels']));
    p=histc(labels,unique(labels))/length(labels);
    effectiven(datai,3,iter) = exp(-sum(p.*log(p)));
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    
    macf1(datai,3,iter)=table2array(f1s(1,1));
    micf1(datai,3,iter)=table2array(f1s(1,2));
    numtables(datai,3,iter) = length(unique(labels));
    slikelihood(:,iter) = readMyMat([data '.likelihood']);
    
    subplot(2,2,3);
    cla;
    scatter(X(:,1),X(:,2),40,labels,'.')
    title([ 'IGMM Slice Sampler: ' num2str(macf1(datai,3,iter))]);
    
%         for j=1:(max(llabels)+1)
%         if (table(j).npoints > 2)
%             sigma = table(j).cholsigma'*table(j).cholsigma;
%             plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0.5 0.5 0.5],2,0.5);
%             %plot_gaussian_ellipsoid(table(j).mu(1:2),sigma(1:2,1:2),'-',[0 0 0],1,1);
%         end
%         end
    
% 
%     tic;
%     [labels,E]=run_dpgmm_subclusters(X', 10, false, 8, false, false, 1, 500, 500);
%     labels = align_labels(labels);
%     p=histc(labels,unique(labels))/length(labels);
%     effectiven(datai,2,iter) = exp(-sum(p.*log(p)));
%     %labels=run_dpgmm_fsd(X',1,false,8,1,40,40);
%     elapsed_time(datai,2,iter)=toc;
%     f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
%     macf1(datai,2,iter)=table2array(f1s(1,1));
%     micf1(datai,2,iter)=table2array(f1s(1,2));
%     numtables(datai,2,iter) = length(unique(labels));
%     changlikelihood(:,iter)=E;
%     subplot(2,2,2);
%     scatter(X(:,1),X(:,2),40,labels,'.')
%     title([ 'JChang Sampler: ' num2str(macf1(datai,2,iter))]);
%     
    
    
    burn_in = '300';
    step = '20';
    dpm_createBinaryFiles(prefix,X,Psi,mu0,m,k0,gam,1);
    cmd = ['dpm64.exe ',data,' ',prior,' ',params,' ',num_sweeps,' ', burn_in,' ',prefix,' ',step];
    fprintf(1,'\nDPGMM is running...\n');
    tic;
    system(cmd);
    elapsed_time(datai,4,iter)=toc;

    [tables customers klabels]=dpm_readOutput(prefix);
    labels = align_labels(klabels);
    f1s=evaluationTable(Y(Y~=0),labels(Y~=0));
    macf1(datai,4,iter)=table2array(f1s(1,1));
    micf1(datai,4,iter)=table2array(f1s(1,2));
    numtables(datai,4,iter) = length(unique(labels));
    p=histc(labels,unique(labels))/length(labels);
    effectiven(datai,4,iter) = exp(-sum(p.*log(p)));
    
    subplot(2,2,4);
    scatter(X(:,1),X(:,2),40,labels,'.')
    title([ 'Collapsed Sampler: ' num2str(macf1(datai,4,iter))]);
    
    
    %subplot(2,2,4);
    %plot([likelihood(:,iter) slikelihood(:,iter) changlikelihood(1:size(likelihood,1),iter)],'linewidth',2);
    %title('Likelihoods');
    %legend(['hetero';'slice ';'jchang'],'Location','southeast');
    %colormap('lines');
    %print(cell2mat(strcat(prefix,'\plots\',names(datai),'_',num2str(iter))),'-depsc');
    end
    subplot(1,1,1);
    

    plot([mean(likelihood,2)],'k:','linewidth',3);hold on;
    plot([mean(slikelihood,2)],'linewidth',3);hold on;
    %plot([mean(changlikelihood(1:size(likelihood,1),:),2)],'r--','linewidth',3);hold off;

    %title('Likelihoods','FontSize',18);
    h=legend(['PPCG';'PSS ';'SUBC'],'Location','southeast');
    set(h,'FontSize',24);
    set(gca,'FontSize',18);
    colormap('lines');
    ylabel('Log Likelihood')
    xlabel('iteration')
    print(cell2mat(strcat(prefix,'\plots\',names(datai),'_',num2str(iter))),'-depsc');
end