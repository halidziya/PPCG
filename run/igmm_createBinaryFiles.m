function igmm_createBinaryFiles(filename,X,Psi,mu0,m,kappa,gamma) % No labels
    if ~exist('data','dir')
        mkdir('data');
    end
    writeMat([ filename '.matrix'],X,'double');
    writeMat([ filename '_params.matrix'],[size(mu0,2) m kappa gamma],'double'); % Aspire format
    writeMat([ filename '_prior.matrix'],[Psi;mu0],'double');
end


function writeMat(filename,mat,prec)
[n m]=size(mat);
file=fopen(filename,'w');
fwrite(file,n,'int');
fwrite(file,m,'int');
fwrite(file,mat',prec);
fclose(file);
end