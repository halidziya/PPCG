function [files,names]=igmm_datasets(dir)

datafolders=ls(dir);
names={};
files={};
for datai=1:size(datafolders,1)
   filename = [dir '\' strtrim(datafolders(datai,:)) '\readData.m'];
   if (exist(filename,'file'));
       files(length(files)+1)={filename};
       names(length(files)) = {strtrim(datafolders(datai,:))};
   end
end
end