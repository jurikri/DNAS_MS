%=====================================================
%Coplyleft, Jan. 28 2019
%Lee Suk-Ho M.D. Ph.D.
%Department of Physiology, Seoul National Univ. College of Medicine
%=====================================================

clear;
%=====================================================
%Given Parameters
%=====================================================
Vrest=-63;
E_K = -80;
E_Na = 50;
E_L = 0;
Cm = 25;
G_K = 28;
G_Na = 150;
G_L = 0.1;

Time=70; %ms
deltaT=0.01;
tw=(0:deltaT:Time)';
Ts1 = 10; Tdur = 5;  Ts2 = Ts1+Tdur;
Ts1 = floor(Ts1/deltaT); Ts2 = floor(Ts2/deltaT);
nt = numel(tw);

Istim=zeros(nt,1);
Istim(Ts1:Ts2,1)=100; % % micro amp current

%=====================================================
%Set Initial values
%=====================================================
dV = 0; 

alpha_n=0.01*(10-dV)/(exp((10-dV)/10)-1);
beta_n=0.125*exp(-dV/80);

alpha_m=0.1*(25-dV)/(exp((25-dV)/10)-1);
beta_m=4*exp(-dV/18);

alpha_h=0.07*exp(-dV/20);
beta_h=1/(exp((30-dV)/10)+1);

ninf=alpha_n/(alpha_n+beta_n);
minf=alpha_m/(alpha_m+beta_m);
hinf=alpha_h/(alpha_h+beta_h);

%====================================================
%Defining arena for simulation
%====================================================
Vm = Vrest*ones(nt, 1);
n = ninf*ones(nt, 1);
m = minf*ones(nt, 1);
h = hinf*ones(nt, 1);

%====================================================
% Simulation Main
%====================================================
for i=1:numel(tw)-1

    dV = Vm(i) - Vrest;
    alpha_n(i) = .01 * ( (10-dV) / (exp((10-dV)/10)-1) );
    beta_n(i) = .125*exp(-dV/80);
    alpha_m(i) = .1*( (25-dV) / (exp((25-dV)/10)-1) );
    beta_m(i) = 4*exp(-dV/18);
    alpha_h(i) = .07*exp(-dV/20);
    beta_h(i) = 1/(exp((30-dV)/10)+1);

    %calc. currents
    I_Na = G_Na * (m(i)^3) *h(i) * (Vm(i)-E_Na);
    I_K = G_K * (n(i)^4) * (Vm(i)-E_K);
    I_L = G_L*(Vm(i)-E_L);
    Itot = Istim(i) - I_K - I_Na - I_L ;

    %calc. derivatives
    Vm(i+1) = Vm(i) + deltaT*Itot/Cm;
    n(i+1) = n(i) + deltaT*(alpha_n(i) *(1-n(i)) - beta_n(i) * n(i));
    m(i+1) = m(i) + deltaT*(alpha_m(i) *(1-m(i)) - beta_m(i) * m(i));
    h(i+1) = h(i) + deltaT*(alpha_h(i) *(1-h(i)) - beta_h(i) * h(i));

end

%===================================================
% Plot Vm
%===================================================
figure(1); clf;
plot(tw,Vm,'LineWidth',2)
axis([0 Time -70 30]);
hold on
legend({'Vm'})
ylabel('mV')
xlabel('ms')
title('Vm');

%====================================================
%Plot Gate dynamics 
%====================================================
figure(2); clf;
p1 = plot(tw, n,'b','LineWidth',2);
hold on
p2 = plot(tw, m,'m','LineWidth',2);
hold on
p3 = plot(tw,h,'g','LineWidth',2);
legend([p1, p2, p3], 'n', 'm', 'h')
ylabel('gate')
xlabel('time (ms)')
title('gate dynamics')