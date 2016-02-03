function [ tables customers labels]=dpm_readOutput(results_dir)
if (~exist('results_dir','var')) %Default value
    results_dir = ''; 
end

tablesfile=fopen([results_dir 'tables.bin']);
customersfile = fopen([results_dir 'customers.bin']);
labelsfile = fopen([results_dir 'Labels.matrix']);

ntable=fread(tablesfile,1,'int');
fprintf(1,'Tables %d',ntable);
tables = readTables(tablesfile,ntable);

ncust=fread(customersfile,1,'int');
customers = readCust(customersfile,ncust);

labels = readMat(labelsfile);

fclose(tablesfile);
fclose(customersfile);
fclose(labelsfile);
end


function mat=readMat(file)
    r = fread(file,1,'int');
    d = fread(file,1,'int');
    mat=fread(file,r*d,'double');
    mat = reshape(mat,d,r); % Row major , column major difference
    % triangle = fread(file,1,'int');
end

function stut=readStut(file)
    stut.eta = fread(file,1,'double');
    stut.normalizer = fread(file,1,'double');
    stut.coef1 = fread(file,1,'double');
    stut.mu = readMat(file);
    stut.cholSigma = readMat(file);
end

function custs=readCust(file,n)
    cust.likelihood0 = fread(file,1,'double');
    cust.tableid = fread(file,1,'int');
    cust.data = readMat(file);
    
    custs =  repmat(cust,1,n); % Just to speed up
    
    for i=2:n
        cust.likelihood0 = fread(file,1,'double');
        cust.tableid = fread(file,1,'int');
        cust.data = readMat(file);
        custs(i) = cust;
    end
end

function tables=readTables(file,n)
  
    table.tableid = fread(file,1,'int');
    table.npoints = fread(file,1,'int');
    table.likelihood = fread(file,1,'double');
    table.scatter = readMat(file);
    table.mean = readMat(file);
    table.dist = readStut(file);
    
    
    tables = repmat(table,1,n); % Speed up
    
    for i=2:n
        table.tableid = fread(file,1,'int');
        table.npoints = fread(file,1,'int');
        table.likelihood = fread(file,1,'double');
        table.scatter = readMat(file);
        table.mean = readMat(file);
        table.dist = readStut(file);
        tables(i) = table;
    end
    
end