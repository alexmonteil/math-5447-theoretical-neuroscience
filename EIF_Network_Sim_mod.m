%%
% This code runs the LIF model with a periodic input
% Designed to explore n:1 entrainment 
%
% Here, V(t) is normalized to be between 0 and 1. If V(t) >= 1, it omits a
% spike and is reset to 0.
%%
clear; close all; clc;

% setting the seed for reproducibility
rng(123);

%%
% Model parameters

% Number of external, exc and inh neurons
Nx=800;
Ne=800;
Ni=200;

% Connection strengths
jex=30.0;
jix=20.0;
jee=12.0;
jei=-45.0;
jie=35.0;
jii=-70.0;

% Connection probabilities
pex=.2;
pix=.2;
pee=.2;
pei=.2;
pie=.2; 
pii=.2;

% External neuron firing rates
%rix
rix=5/1000;
rex=5/1000;

% EIF neuron parameters
taum=10; 
EL=-72;
Vre=-72; 
VT=-55;
D=2;
Vth=0;

% Synaptic time constants
taux=8;
taue=6;
taui=4;

% Discretized time
T=200;
dt=.1;
time=0:dt:T;
Nt=length(time);

%% Input drive and the connectivity matrix
% Generate external spike trains as Poisson processes
Six=binornd(1,rix*dt,Nx,length(time))/dt;
Sex=binornd(1,rex*dt,Nx,length(time))/dt;

% Connection matrices
Jex=jex*binornd(1,pex,Ne,Nx);
Jee=jee*binornd(1,pee,Ne,Ne);
Jei=jei*binornd(1,pei,Ne,Ni);
Jix=jix*binornd(1,pix,Ni,Nx);
Jie=jie*binornd(1,pie,Ni,Ne);
Jii=jii*binornd(1,pii,Ni,Ni);

%%
[Se, Si] = EIF_Network_Euler_mod(Ne, Ni,VT,EL,time,dt,D,taum,taux, taue, taui, ...
    Vth, Vre,Sex, Six,Jex,Jee,Jei,Jix,Jie,Jii);

%% Plot the results

% Get spike times and neuron indices
[SxIndices,SxTimes]=find(Six>0); 
SxTimes=SxTimes*dt;
[SeIndices,SeTimes]=find(Se>0);
SeTimes=SeTimes*dt;
[SiIndices,SiTimes]=find(Si>0);
SiTimes=SiTimes*dt;


% Compute population-averaged rates across time
re=mean(Se);
ri=mean(Si);

% Compute excitatory-network average rate over time
avg_re = mean(re(1001:end));

% Compute inhibitory-network average rate over time
avg_ri = mean(ri(1001:end));

% Compute network average rate across both excitatory and inhibitory over
% time
network_avg_rate = (avg_re * Ne + avg_ri * Ni) / (Ne + Ni);

% Display results
fprintf('Average E-rate: %.2f Hz\n', avg_re * 1000);
fprintf('Average I-rate: %.2f Hz\n', avg_ri * 1000);
fprintf('Overall network average rate: %.2f Hz\n', network_avg_rate * 1000);

% Smooth the population-averaged rates
sigma=6;
s=(-3*sigma:dt:3*sigma);
k=exp(-(s.^2)/(2*sigma^2)); % Gaussian kernel
k(s<0)=0; % Make it causal
k=k/(sum(k)*dt);
re=conv(re,k,'same')*dt;
ri=conv(ri,k,'same')*dt;

%size(re)
disp(mean(re(1001:end)))
disp(length(Se))
%disp(re)
disp(mean(re))
%disp(ri)

%%
subplot(2,3,1)
plot(SxTimes,SxIndices,'k.','markersize',10)
subplot(2,3,2)
plot(SeTimes,SeIndices,'.','markersize',10,'color',[0.8660 0.3290 0])
subplot(2,3,5)
plot(time,re*1000,'color',[0.8660 0.3290 0])


subplot(2,3,3)
plot(SiTimes,SiIndices,'.','markersize',10,'color',[0.066 0.443 0.7450])

subplot(2,3,6)
plot(time,ri*1000,'color',[0.066 0.443 0.7450])
