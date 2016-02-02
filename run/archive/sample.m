% load('ToyData_I2GMM_journal_3rare_3normal_3D_final.mat')
% load fisheriris
% X=meas;
% ul=unique(species);
% Y=zeros(150,1);
% for i=1:length(species)
%     Y(i) = find(strcmp(ul,species(i)));
% end

%load ionosphere
%Y=strcmp(Y,'g')+1;

X=[ mvnrnd([0,0],eye(2)/2,300); mvnrnd([1,2],[0.6 -0.4;-0.4 0.6],300); mvnrnd([1,-2],[0.6 0.4;0.4 0.6],300)];
Y=[ones(300,1);2*ones(300,1);3*ones(300,1)];

[labels M]=UncollapsedSampler(X,Y);