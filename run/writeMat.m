function writeMat(filename,mat,prec)
[n m]=size(mat);
file=fopen(filename,'w');
fwrite(file,n,'int');
fwrite(file,m,'int');
fwrite(file,mat',prec);
fclose(file);
end

