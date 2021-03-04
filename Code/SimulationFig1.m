
%Initialization
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

K3t = [0 0.05 0.1 0.2]; %Range of Error Vector Magnitudes for transmission at the relay
K3r = K3t; %Range of Error Vector Magnitudes for reception at the relay (same as for transmission)


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

PoutSimulation1 = zeros(length(P1),length(K3t)); %Outage probability at transmitter 1: Based on Monte-Carlo simulations
PoutSimulation2 = zeros(length(P1),length(K3t)); %Outage probability at transmitter 2: Based on Monte-Carlo simulations

PoutAnalytic1 = zeros(length(P1),length(K3t)); %Outage probability at transmitter 1: Based on analytic expressions
PoutAnalytic2 = zeros(length(P1),length(K3t)); %Outage probability at transmitter 1: Based on analytic expressions

PoutAsymptoticLimit1 = zeros(1,length(K3t)); %Outage probability at transmitter 1: Asymptotic limit as transmit powers go to infinity
PoutAsymptoticLimit2 = zeros(1,length(K3t)); %Outage probability at transmitter 2: Asymptotic limit as transmit powers go to infinity


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
        
        %Compute the average number of outages by Monte-Carlo simulations
        PoutSimulation1(p,k) = sum(SNDR1(p,:,k) < x) /nbrOfRealizations;
        PoutSimulation2(p,k) = sum(SNDR2(p,:,k) < x) /nbrOfRealizations;
        
        
        %Part 2: Analytic results
        
        %Compute the different terms that show up in Eq. (10)
        expArgument1 = - ( x/(1-c*x) * (a1/Omega2+b1/Omega1) + x*(1+c*x)/(1-c*x)^2 *b1/Omega2 *P1(p)/P2(p)); %Argument of first exponential function in Eq. (10), for transmitter 1
        expArgument2 = - ( x/(1-c*x) * (a2/Omega1+b2/Omega2) + x*(1+c*x)/(1-c*x)^2 *b2/Omega1 *P2(p)/P1(p)); %Argument of first exponential function in Eq. (10), for transmitter 2
        numerator1 = (x+x^2)/(1-c*x)^2*N1*N3/(Omega1*Omega2*P2(p)*P3(p)) + x^2/(1-c*x)^3 * b1^2*P1(p)/(Omega1*Omega2*P2(p)); %Numerator of the expression in the square roots in Eq. (10), for transmitter 1
        numerator2 = (x+x^2)/(1-c*x)^2*N2*N3/(Omega1*Omega2*P1(p)*P3(p)) + x^2/(1-c*x)^3 * b2^2*P2(p)/(Omega1*Omega2*P1(p)); %Numerator of the expression in the square roots in Eq. (10), for transmitter 2
        denominator1 = 1 + c*x/(1-c*x) * P1(p)*Omega1/(P2(p)*Omega2); %Denominator of the expression in the square roots in Eq. (10), for transmitter 1
        denominator2 = 1 + c*x/(1-c*x) * P2(p)*Omega2/(P1(p)*Omega1); %Denominator of the expression in the square roots in Eq. (10), for transmitter 2
        
        %Compute outage probability analytically using Proposition 1
        if c<1/x
            PoutAnalytic1(p,k) = 1 - exp(expArgument1) *2*sqrt(numerator1/denominator1) * besselk(1,2*sqrt(numerator1*denominator1));
            PoutAnalytic2(p,k) = 1 - exp(expArgument2) *2*sqrt(numerator2/denominator2) * besselk(1,2*sqrt(numerator2*denominator2));
        else
            PoutAnalytic1(p,k) = 1;
            PoutAnalytic2(p,k) = 1;
        end
        
    end
    
end


%Compute the asymptotic outage probability according to Corollary 1
for k = 1:length(K3t)
    c = K3t(k)^2 + K3r(k)^2 + K3t(k)^2 *K3r(k)^2; %Parameter c
    
    %Compute asymptotic limits according to Eq. (8)
    if c < 1/x
        PoutAsymptoticLimit1(k) = Omega1*c*x/(Omega2+c*x*(Omega1-Omega2));
        PoutAsymptoticLimit2(k) = Omega2*c*x/(Omega1+c*x*(Omega2-Omega1));
    else
        PoutAsymptoticLimit1(k) = 1;
        PoutAsymptoticLimit2(k) = 1;
    end
end


%Plot outage probability results for transmitter 1
figure(1); hold on; box on;

for k = 1:length(K3t)
    plot(P1_dB,PoutSimulation1(:,k),'k*','Markersize',6); %Markers based on Monte Carlo simulations
    plot(P1_dB,PoutAnalytic1(:,k),'k','Linewidth',1); %Lines based on analytic results
    plot([P1_dB(1) P1_dB(end)],PoutAsymptoticLimit1(k)*ones(1,2),'--k','Linewidth',1.5); %Asymptotic limit
    text(P1_dB(end-3),PoutAnalytic1(end-3,k)+0.05,['\kappa = ' num2str(K3t(k))]);
end
axis([P1_dB(1) P1_dB(end) 0 1.1]);

legend('Simulation', 'Analytical (Exact)','Analytical (Asympotic)','Location','Best');
%xlabel('Transmit Powers P_1 and P_2 [dB]');
%ylabel('Outage Probability at Transmitter Node 1');
axis([0 50 -0.02 1.1])



%Plot outage probability results for transmitter 2.
%This figure is not shown in the paper since it is identical to the
%results for transmitter 1
figure(2); hold on; box on;

for k = 1:length(K3t)
    plot(P1_dB,PoutSimulation2(:,k),'k*','Markersize',6); %Markers based on Monte Carlo simulations
    plot(P1_dB,PoutAnalytic2(:,k),'k','Linewidth',1); %Lines based on analytic results
    plot([P1_dB(1) P1_dB(end)],PoutAsymptoticLimit2(k)*ones(1,2),'--k','Linewidth',1.5); %Asymptotic limit
    text(P1_dB(end-3),PoutAnalytic2(end-3,k)+0.05,['\kappa = ' num2str(K3t(k))]);
end
axis([P1_dB(1) P1_dB(end) 0 1.1]);
title('Secrecy outage probability against interference power for case 1','FontWeight','Normal');

legend('Simulation', 'Analytical (Exact)','Analytical (Asympotic)','Location','Best');
xlabel('Interference Power (dB) ');
ylabel('Secrecy Outage Probability');
axis([0 50 -0.02 1.1])
