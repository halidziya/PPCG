function [ tables,labels]=igmm_readOutput(filename)
file=fopen(filename);
ntables=fread(file,1,'int');
tables = [];
for i=1:ntables
    table=readTable(file);
    tables = [tables table];
end
labels = readMat(file);
fclose(file);
end

function mat=readMat(file)
    r = fread(file,1,'int');
    d = fread(file,1,'int');
    mat=fread(file,r*d,'double');
    mat = reshape(mat,d,r);
    % triangle = fread(file,1,'int');
end

function table=readTable(file)
        table.npoints = fread(file,1,'int');
        table.mu = readMat(file);
        table.cholsigma = readMat(file);
        table.scatter = readMat(file);
        table.sum = readMat(file);
end
