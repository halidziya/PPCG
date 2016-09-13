function slice_createBinaryFiles(filename,X,Psi,mu0,m,kappa,gamma) % No labels
    writeMat([ filename '.matrix'],X,'double');
    writeMat([ filename '_params.matrix'],[size(mu0,2) m kappa gamma]','double'); % Column vector
    writeMat([ filename '_psi.matrix'],Psi,'double');
    writeMat([ filename '_mean.matrix'],mu0,'double');
end


