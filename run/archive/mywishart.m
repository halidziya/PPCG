function w=mywishart(d,df)
   v=zeros(d,1);
   n=normrnd(0,1,d,d);
   w=zeros(d,d);
   for i=1:d
        v(i) = chi2rnd(df-i+1);
        w(i,i)=v(i)+sum(n(i,1:(i-1)).^2);
        for j=1:d
            if (i<j)
            w(i,j) = n(i,j)*sqrt(v(i))+n(i,1:(i-1))*n(1:(i-1),j);             
            end
        end  
   end
   
   for i=1:d
       for j=1:(i-1)
           w(i,j) = w(j,i);
       end
   end
end