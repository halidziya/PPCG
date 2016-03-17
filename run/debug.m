
a=fopen('C:\Users\hzyereba\Desktop\UncollapsedSampler\Qtigmm\tlikelihood.matrix','rb')
tlikelihood = readMyMat(a);
fclose(a);

a=fopen('C:\Users\hzyereba\Desktop\UncollapsedSampler\Qtigmm\lastlikelihood.matrix','rb')
likelihood = readMyMat(a);
fclose(a);
scatter(X(:,1),X(:,2),exp(tlikelihood-likelihood)*150)

