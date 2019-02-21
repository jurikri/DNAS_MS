%=====================================================
%Coplyleft, Jan. 28 2019
%Lee Suk-Ho M.D. Ph.D.
%Department of Physiology, Seoul National Univ. College of Medicine
%=====================================================

clear; close all;
% Load the training pattern
img = imread('s', 'bmp');
Npc = numel(img); % # of PCs = 3600 # pix 갯수 

% Constants
time = 30; % in ms
c = 0.25; % connectivity btwn PCs
g1 = 0.5; % network activity-dependence of threshold 
b1 = 0.1; % partial cue to the original one # 전체 s 중에서, partial cue를 줄 때, 10%만 주겠다. 

%Izhs param
vrest = -65; vret = -58;
v = vrest*ones(Npc,1); 
aiz=0.02;  biz = 0.25;  ud = 4; 
WtHigh = 2;
WtLow = WtHigh/2;

%% First training

img = imread('s', 'bmp');
Npc = numel(img); % # of PCs = 3600 # pix 갯수 
img = min(1,-1*(double(img)-255));

W = rand(Npc,Npc); % W = structural connection matrix
W = ceil(W + (c-1));

sum(sum((W == 0)))/(size(W,1)*size(W,2))

Z=zeros(Npc,1);  %training pattern
Z(:,1)=img(:);

% Jh = Hebbian connectivity
Jh = Z*Z';
Jh = min(WtLow,Jh);

%Disable autaptic connections
for idx=1:Npc
    Jh(idx,idx) = 0;
end

% general = 1; % ??? ?? ??? ?? ?? ?


%% arc function
degradation = 1;
tSET = 2;

%%
lr_list = find(Z(:,1)==1);
not_lr_list = find(Z(:,1)==0);
WJ = max(WtLow*W, WtHigh*W.*Jh); % Jh = Hebbian connectivity
for X = 1:size(lr_list,1) % self: X
    for Y = 1:size(not_lr_list,1)
        WJ(lr_list(X,1),not_lr_list(Y,1)) = WJ(lr_list(X,1),not_lr_list(Y,1))*degradation;
    end
end

figure()
imagesc(WJ)

sum(sum((WJ == 0)))/(size(WJ,1)*size(WJ,2))

%Syn Wt matrix
%WJ = SynWt*W.*Jh;



%% Second training: t import
img2 = imread('t', 'bmp');
img2 = min(1,-1*(double(img2)-255));

% W = rand(Npc,Npc);
% W = ceil(W + (c-1));

Z2=zeros(Npc,1);  %training pattern
Z2(:,1)=img2(:);
sum(sum(img))
sum(sum(img2))
% Jh = Hebbian connectivity
Jh2 = Z2*Z2';
Jh2 = min(WtLow,Jh2);
(sum(sum(Jh2==0))+sum(sum(Jh2==1)))/3600
%Disable autaptic connections
for idx=1:Npc
    Jh2(idx,idx) = 0;
end
(sum(sum(Jh2==0))+sum(sum(Jh2==1)))/3600
% Jh = (Jh + Jh2); % key

% Jh_tmp = Jh;
% figure(1)
% imagesc(WJ)

% max ???? ???? ??? ??? ????? ??, ?? ?????? ???? ??? ? 

lr_list2 = find(Z2(:,1)==1);
not_lr_list2 = find(Z2(:,1)==0);

WJ2 = WtHigh*(W.*Jh2); % Jh = Hebbian connectivity

sum(sum(abs(WJ2(W == 0))))
sum(sum(abs(WJ(W == 0))))

sum(sum((W == 0)))/(size(W,1)*size(W,2))
sum(sum((WJ2 == 0)))/(size(WJ2,1)*size(WJ2,2))
sum(sum((Jh2 == 0)))/(size(Jh2,1)*size(Jh2,2))

for row = 1:size(WJ,1)
    for col = 1:size(WJ,2)
        if W(row,col) == 0
            WJ(row,col) = 0 +  WJ2(row,col);
            if WJ2(row,col) ~= 0
                disp('t')
                disp(WJ2(row,col))
            end
        elseif W(row,col) == 1
            WJ(row,col) = max(WJ(row,col),WJ2(row,col));
        elseif WJ(row,col) == 2
            WJ(row,col) = WJ(row,col) + WJ2(row,col);
        else
            WJ(row,col)
        end
    end
end

sum(sum((WJ == 0)))/(size(WJ,1)*size(WJ,2))

% WJ = WJ + WJ2;
lr_list2 = find(Z2(:,1)==1);
not_lr_list2 = find(Z2(:,1)==0);
for X = 1:size(lr_list2,1) % self: X
    for Y = 1:size(not_lr_list2,1)
            WJ(lr_list2(X,1),not_lr_list2(Y,1)) = WJ(lr_list2(X,1),not_lr_list2(Y,1))*degradation;
    end
end
figure(1)
imagesc(WJ)


% figure(2)
figure(2)
imagesc(WJ)

%% Third training: t import
% img3 = (img + img2) >0;
% W = rand(Npc,Npc);
% W = ceil(W + (c-1));
% Z2=zeros(Npc,1);  %training pattern
% Z2(:,1)=img3(:);
% 
% % Jh = Hebbian connectivity
% Jh = Z2*Z2';
% Jh = min(WtLow,Jh);
% 
% %Disable autaptic connections
% for idx=1:Npc
%     Jh(idx,idx) = 0;
% end
% 
% figure()
% imagesc(Jh_tmp)
% 
% figure()
% imagesc(Jh)
% 
% sum(sum((Jh2 - Jh) ~= 0))
% 
% figure(1)
% imagesc(WJ)
% WJ = max(general*W, WtHigh*W.*Jh); % Jh = Hebbian connectivity
% figure(2)
% imagesc(WJ)

%% hebbian connection? ???? High? ??, ??? low? ???.

Isyn = zeros(Npc,time); % h(i,t): exc synaptic inputs to i-th PC at t 
Inet = zeros(Npc, time);
Iinh = zeros(time,1); % Inh. synaptic inputs at t

%Make initial cue
if tSET == 1
    X0 = double(img(:)).*rand(Npc, 1);
%% Make initial seconde cue 
elseif tSET == 2
    X0 = double(img2(:)).*rand(Npc, 1);
end
% 
%%
% figure(1)
% plot(X0)
% tmp1 = X0+(b1-1);
% figure(2)
% plot(tmp1)
X0 = max(0,ceil(X0+(b1-1))); %initial cue % initial cue? ??? 10%? ?????.
% figure(3)
% plot(X0) 

X = zeros(Npc,time); % X(i,t) = 1, if ith PC fired at t.
X(:,1) = X0;

% plot initial partial cue
fx0 = find(X0(:)>0);
nrow = size(img,1); ncol = size(img,2); np = ncol*nrow;
figure(5); clf; axis tight equal  
col = ceil(fx0/nrow);    row = mod(fx0, nrow);    plot(col, row, '.');  axis([0 60 0 60]); 
pause(2); clear fx0;

v(:,1) = vrest; uiz = biz.*v; % potential ??, biz? %Izhs param
Iinh(:,1) = 0; Inet(:,1) = 0; Isyn(:,1) = 0;

%%

% Simulation Main
for t = 1:time-1
    t
	Isyn(:,t) = WJ*X(:,t);  % Isyn(i,t): exc. synaptic inputs to i-th PC at t % X? firing ??,
	Iinh(t) = g1*sum(X(:,t));  % Inh. synaptic inputs
    Inet(:,t) = max(0, Isyn(:,t) - Iinh(t));
%     Iinh(t)
        
    fired=find(v>=30);   
    v(fired)=vret;       % reset of fired
	uiz(fired)=uiz(fired)+ud; 
    X(fired,t+1) = 1;
    
    v=v+0.5*((0.04*v+5).*v+140-uiz+ Inet(:,t));
	v=v+0.5*((0.04*v+5).*v+140-uiz+ Inet(:,t));
    
    uiz=uiz+aiz.*(biz*v-uiz);
    
    hold on
    if mod(t,2)==1       
       pause(0.2)
    end
    col = ceil(fired/nrow);    row = mod(fired, nrow);    plot(col, row, '.');   axis([0 60 0 60]);      
end

























