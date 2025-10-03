function [Se, Si] = EIF_Network_Euler(Ne,Ni,VT,EL,time,dt,D, taum, taux, ...
    taue, taui, Vth, Vre, Sx,Jex,Jee,Jei,Jix,Jie,Jii)


% Initialize membrane potentials randomly
Ve=rand(Ne,1)*(VT-EL)+EL;
Vi=rand(Ni,1)*(VT-EL)+EL;
% Initialize all else as zeros
Iex=zeros(Ne,1);
Iee=zeros(Ne,1);
Iei=zeros(Ne,1);
Iix=zeros(Ni,1);
Iie=zeros(Ni,1);
Iii=zeros(Ni,1);
Se=zeros(Ne,length(time));
Si=zeros(Ni,length(time));

for tt = 1:length(time)-1
    % Euler step for membrane potentials
    Ve=Ve+dt*(-(Ve-EL)+D*exp((Ve-VT)/D)+Iex+Iee+Iei)/taum; 
    Vi=Vi+dt*(-(Vi-EL)+D*exp((Vi-VT)/D)+Iix+Iie+Iii)/taum; 
    
    % Euler step for synaptic currents
    Iex=Iex+dt*(-Iex+Jex*Sx(:,tt))/taux;      
    Iee=Iee+dt*(-Iee+Jee*Se(:,tt))/taue;
    Iei=Iei+dt*(-Iei+Jei*Si(:,tt))/taui;   
    Iix=Iix+dt*(-Iix+Jix*Sx(:,tt))/taux;    
    Iie=Iie+dt*(-Iie+Jie*Se(:,tt))/taue;
    Iii=Iii+dt*(-Iii+Jii*Si(:,tt))/taui;
    
        
    % Find which excitatory neurons spiked.
    Inds= Ve>=Vth;
    % Reset membrane potentials
    Ve(Inds)=Vre;    
    % Store spikes as delta functions
    Se(Inds,tt+1)=1/dt;
        
    % Now do the same for inhibitory neurons
    Inds= Vi>=Vth;
    % Reset membrane potentials
    Vi(Inds)=Vre;
    % Store spikes as delta functions
    Si(Inds,tt+1)=1/dt;
end