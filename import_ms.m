ForCheck = input('시작하려면 아무 숫자를 입력하세요 '); % 개발중에 실수로 시작하여 변수가 날아가는것을 방지하기 위한 체크 구문
clear; clc;

path = 'C:\Users\msbak\Desktop\DANS_team\DANS\NeuralData\project';
%%
raw_A = file_import_spike([path 'A\']);
raw_B = file_import_spike([path 'B\']);



%%
