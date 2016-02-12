function W=myinvwish(sigma,df)
   d=size(sigma,2);
   w=zeros(d,d);
   for i=1:d
        w(i,i)=sqrt(chi2rnd(df-i+1));
        for j=1:d
            if (i<j)
            w(i,j) =normrnd(0,1);             
            end
        end  
   end
   T = chol(sigma)'*inv(w);
   W = T*T';
end