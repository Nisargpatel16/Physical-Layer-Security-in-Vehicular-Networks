
close all;
clear all;


%Two-way relaying with transceiver hardware impairments

R = 5;       %Transmit information rate (per 2 channel uses)
x = 2.^R-1;  %Corresponding outage thresholds

Omega1 = 1;  %Average channel gain between transmitter 1 and relay
Omega2 = 1;  %Average channel gain between transmitter 2 and relay

N1 = 1;      %Normalized noise variance at transmitter 1
N2 = 1;      %Normalized noise variance at transmitter 2
N3 = 1;      %Normalized noise variance at relay

P1_dB = 0:2.5:50;    %SNR between transmitter 1 and relay
P2_dB = P1_dB;       %SNR between transmitter 2 and relay
P1 = 10.^(P1_dB/10); %Corresponding transmit power range at transmitter 1
P2 = 10.^(P2_dB/10); %Corresponding transmit power range at transmitter 2
P3 = 0.5*P2;         %Relay transmit power (half of transmitters)

K3t = [0 0.05 0.1]; %Range of Error Vector Magnitudes for transmission at the relay
K3r = 0.2 - K3t; %Range of Error Vector Magnitudes for reception at the relay (sum of transmission and reception is 0.2)


%Number of realizations in Monte-Carlo simulations
nbrOfRealizations = 1000000;

%Generate channel fading realizations. RHO1 and RHO2 are the squared norms.
h1 = sqrt(Omega1/2)*(randn(1,nbrOfRealizations) + 1i*randn(1,nbrOfRealizations)); %Rayleigh fading for channel from transmitter 1 to relay
h2 = sqrt(Omega2/2)*(randn(1,nbrOfRealizations) + 1i*randn(1,nbrOfRealizations)); %Rayleigh fading for channel from transmitter 2 to relay
RHO1 = abs(h1).^2;  %Fading power of h1
RHO2 = abs(h2).^2;  %Fading power of h2

%Modulation parameters for symbol error rate computation (using Eq. (15))
alpha = 1; %BPSK modulation
beta = 1;  %BPSK modulation


%Placeholders for simulation results
SNDR1 = zeros(length(P1),nbrOfRealizations,length(K3t)); %Signal-to-noise-and-distortion ratio (SNDR) at transmitter 1
SNDR2 = zeros(length(P1),nbrOfRealizations,length(K3t)); %Signal-to-noise-and-distortion ratio (SNDR) at transmitter 2

SERSimulation1 = zeros(length(P1),length(K3t)); %Symbol error rate at transmitter 1: Based on Monte-Carlo simulations
SERSimulation2 = zeros(length(P1),length(K3t)); %Symbol error rate at transmitter 2: Based on Monte-Carlo simulations

SERAnalytic1 = zeros(length(P1),length(K3t));  %Symbol error rate at transmitter 1: Based on analytic expressions
SERAnalytic2 = zeros(length(P1),length(K3t));  %Symbol error rate at transmitter 2: Based on analytic expressions

SERasymptoticLimit1 = zeros(1,length(K3t)); %Symbol error rate at transmitter 1: Asymptotic limit as transmit powers go to infinity
SERasymptoticLimit2 = zeros(1,length(K3t)); %Symbol error rate at transmitter 2: Asymptotic limit as transmit powers go to infinity


%Go through all cases of transmit power at the transmitters and assume that
%the power on the relay is varied in the same way. There are also different
%cases with different levels of transceiver hardware impairments.
for p = 1:length(P1)
    for k = 1:length(K3t)
        
        %Relaying gain as in Eq. (5)
        G = sqrt( P3(p) ./ ( (RHO1*P1(p)+RHO2*P2(p))*(1+K3r(k)^2) + N3 ) );
        
        %Compute the 5 parameters that are defined in Section II, right
        %after Eq. (8).
        a1 = (N3/P2(p)) * (1+K3t(k)^2);
        a2 = (N3/P1(p)) * (1+K3t(k)^2);
        b1 = (N1/P3(p)) * (1+K3r(k)^2);
        b2 = (N2/P3(p)) * (1+K3r(k)^2);
        c = K3t(k)^2 + K3r(k)^2 + K3t(k)^2 *K3r(k)^2;
        
        %Compute SNDRs for different fading realizations, according to Eq. (8)
        SNDR1(p,:,k) = (RHO1.*RHO2) ./ ( RHO1.^2*(P1(p)/P2(p))*c + RHO1.*RHO2 *c + RHO2*b1 + RHO1*(a1+(P1(p)/P2(p))*b1) + N1*N3/(P2(p)*P3(p)) );
        SNDR2(p,:,k) = (RHO1.*RHO2) ./ ( RHO2.^2*(P2(p)/P1(p))*c + RHO1.*RHO2 *c + RHO1*b2 + RHO2*(a2+(P2(p)/P1(p))*b2) + N2*N3/(P1(p)*P3(p)) );
        
        
        %Part 1: Monte-Carlo simulations
        
        %Compute SERs for different fading realizations, according to Eq. (15)
        SERSimulation1(p,k) = mean(alpha*qfunc(sqrt(2*beta*SNDR1(p,:,k)))); %Simulated SER at transmitter 1
        SERSimulation2(p,k) = mean(alpha*qfunc(sqrt(2*beta*SNDR2(p,:,k)))); %Simulated SER at transmitter 2
        
        
        %Part 2: Analytic results
        
        %Compute outage probability analytically by computing Eq. (16) the
        %integration numerically
        SERAnalytic1(p,k) = quadgk(@(x) functionSERintegrand(x,a1,b1,c,Omega1,Omega2,P1(p),P2(p),P3(p),N1,N3,alpha,beta),0,Inf);
        SERAnalytic2(p,k) = quadgk(@(x) functionSERintegrand(x,a2,b2,c,Omega2,Omega1,P2(p),P1(p),P3(p),N2,N3,alpha,beta),0,Inf);
        
    end
    
end


%Compute the asymptotic SER according to Eq. (17). This formula only
%holds when the average channel gains are the same.
if Omega1 == Omega2
    for k = 1:length(K3t)
        c = K3t(k)^2 + K3r(k)^2 + K3t(k)^2 *K3r(k)^2; %Parameter c
        
        %Limit in Eq. (17). Note: gammainc in Matlab has an
        %unconventional normalization factor, therefore the expression
        %in Eq. (17) has been multiplied by gamma(3/2).
        SERasymptoticLimit1(k) = (alpha*c)/(2*beta*sqrt(pi)) * gammainc(beta/c,3/2) *gamma(3/2) + (alpha/2) * erfc(sqrt(beta/c));
        SERasymptoticLimit2(k) = SERasymptoticLimit1(k);
    end
end


%Plot symbol error rate results for transmitter 1.
figure(1); hold on; box on; grid on;

for k = 1:length(K3t)
    plot(P1_dB,SERSimulation1(:,k),'k*','Markersize',6); %Markers based on Monte Carlo simulations
    plot(P1_dB,SERAnalytic1(:,k),'k','Linewidth',1.0); %Lines based on analytic results
    plot([P1_dB(1) P1_dB(end)],SERasymptoticLimit1(k)*ones(1,2),'--k','Linewidth',1.5); %Asymptotic limit
end

text(P1_dB(end-5),SERAnalytic1(end-5,1)*1.2,['\kappa_3_t = ' num2str(K3t(1)) ', \kappa_3_r = ' num2str(K3r(1))]);
text(P1_dB(end-5),SERAnalytic1(end-5,2)*1.2,['\kappa_3_t = ' num2str(K3t(2)) ', \kappa_3_r = ' num2str(K3r(2))]);
text(P1_dB(end-5),SERAnalytic1(end-5,3)*0.8,['\kappa_3_t = ' num2str(K3t(3)) ', \kappa_3_r = ' num2str(K3r(3))]);

legend('Simulation', 'Numerical (Exact)', 'Analytical (Asympotic)','Location','NorthEast');
xlabel('Transmit Powers P_1 and P_2 [dB]');
ylabel('Symbol Error Rate at Transmitter Node 1');

set(gca,'yscale','log')
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')
set(gca,'XTick',0:10:50)
set(gca,'YTick',[1e-2 1e-1])

axis([0 50 1e-3 1])


%Plot symbol error rate results for transmitter 2.
%This figure is not shown in the paper since it is identical to the
%results for transmitter 1
figure(2); hold on; box on; grid on;

for k = 1:length(K3t)
    plot(P1_dB,SERSimulation2(:,k),'k*','Markersize',6); %Markers based on Monte Carlo simulations
    plot(P1_dB,SERAnalytic2(:,k),'k','Linewidth',1.0); %Lines based on analytic results
    plot([P1_dB(1) P1_dB(end)],SERasymptoticLimit2(k)*ones(1,2),'--k','Linewidth',1.5); %Asymptotic limit
end

text(P1_dB(end-5),SERAnalytic2(end-5,1)*1.2,['\kappa_3_t = ' num2str(K3t(1)) ', \kappa_3_r = ' num2str(K3r(1))]);
text(P1_dB(end-5),SERAnalytic2(end-5,2)*1.2,['\kappa_3_t = ' num2str(K3t(2)) ', \kappa_3_r = ' num2str(K3r(2))]);
text(P1_dB(end-5),SERAnalytic2(end-5,3)*0.8,['\kappa_3_t = ' num2str(K3t(3)) ', \kappa_3_r = ' num2str(K3r(3))]);
title('Secrecy outage probability against interference power for case 2','FontWeight','Normal');

legend('Simulation', 'Numerical (Exact)', 'Analytical (Asympotic)','Location','NorthEast');
xlabel('Interference Power [dB]');
ylabel('Secrecy Outage Probability');

set(gca,'yscale','log')
set(gca,'XTickMode','manual')
set(gca,'YTickMode','manual')
set(gca,'XTick',0:10:50)
set(gca,'YTick',[1e-2 1e-1])

axis([0 50 1e-3 1])