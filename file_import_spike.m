function [msCell, msName] = file_import_spike(path2,s)
ix = 1
file2 = dir([char([path2]) '\*.*'])
for i = 1:size(file2,1)
    try
       tmp = file2(i).name;
       filepath = [path2 tmp];
       msCell{ix} = dlmread(filepath, ',', 13, s);
       msName{ix} = [path2 tmp];
       ix = ix+1;
    catch
    end
end

$ msName에 이름이 저장됩니다
$ s는 추출하고자 하는 column의 시작 위치입니다. 
