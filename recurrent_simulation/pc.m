%=====================================================
%Coplyleft, Jan. 28 2019
%Lee Suk-Ho M.D. Ph.D.
%Department of Physiology, Seoul National Univ. College of Medicine
%=====================================================

clear;
% Load the training pattern
img = imread('s', 'bmp');
Npc = numel(img); % # of PCs = 3600 # pix 갯수 
img = min(1,-1*(double(img)-255));

% Constants
time = 300; % in ms
c = 0.25; % connectivity btwn PCs
g1 = 0.5; % network activity-dependence of threshold 
b1 = 0.1; % partial cue to the original one # 전체 s 중에서, partial cue를 줄 때, 10%만 주겠다. 


%Izhs param
vrest = -65; vret = -58;
v = vrest*ones(Npc,1); 
aiz=0.02;  biz = 0.25;  ud = 4; 
WtHigh = 2;
WtLow = WtHigh/2; 

% W = structural connection matrix
W = rand(Npc,Npc);
W = ceil(W + (c-1));

Z=zeros(Npc,1);  %training pattern
Z(:,1)=img(:);

% Jh = Hebbian connectivity
Jh = Z*Z';
Jh = min(1,Jh);

%Disable autaptic connections
for idx=1:Npc
    Jh(idx,idx) = 0;
end

%Syn Wt matrix
%WJ = SynWt*W.*Jh;
WJ = max(WtLow*W, WtHigh*W.*Jh); % Jh = Hebbian connectivity
% hebbian connection? ???? High? ??, ??? low? ???.

Isyn = zeros(Npc,time); % h(i,t): exc synaptic inputs to i-th PC at t 
Inet = zeros(Npc, time);
Iinh = zeros(time,1); % Inh. synaptic inputs at t

%Make initial cue
X0 = double(img(:)).*rand(Npc, 1);
% figure(1)
% plot(X0)
tmp1 = X0+(b1-1);
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
figure(1); clf; axis tight equal  
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
    v_vis(:,t) = v;
    
    uiz=uiz+aiz.*(biz*v-uiz);
    
    hold on
    if mod(t,2)==1       
       pause(0.2)
    end
    col = ceil(fired/nrow);    row = mod(fired, nrow);    plot(col, row, '.');   axis([0 60 0 60]);      
end

























