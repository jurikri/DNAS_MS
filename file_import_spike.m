function msCell = file_import_spike(path2)
ix = 1
file2 = dir([char([path2]) '\*.*'])
for i = 1:size(file2,1)
    try
       tmp = file2(i).name;
       filepath = [path2 tmp];
       msCell{ix} = dlmread(filepath, ',', 13, 15);
       ix = ix+1;
    catch
    end
end