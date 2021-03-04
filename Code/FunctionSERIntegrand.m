function integrand = functionSERintegrand(x,a_i,b_i,c,Omega_i,Omega_ri,Pi,P_ri,P3,N_i,N3,alpha,beta)
%Compute the integrand of Eq. (16) in the paper:
%
%Michail Matthaiou, Agisilaos Papadogiannis, Emil Björnson, Merouane
%Debbah, "Two-way Relaying under the Presence of Relay Transceiver Hardware
%Impairments," IEEE Communications Letters, vol. 17, no. 6, pp. 1136-1139,
%June 2013.
%
%Download article: http://arxiv.org/pdf/1307.2923
%
%This is version 1.0 (Last edited: 2014-03-20)
%

%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.
%
%
%The input variables follow the notation from Eq. (16) in the paper

%Compute factor in Eq. (16) that is multiplied with the outage probability
factorIntegrand = (alpha*sqrt(beta/pi)/2) * (exp(-beta*x)./sqrt(x));

%Compute the different terms that show up in Eq. (10)
expArgument = - ( x./(1-c*x) * (a_i/Omega_ri+b_i/Omega_i) + x.*(1+c*x)./(1-c*x).^2 * (b_i/Omega_ri *Pi/P_ri) ); %Argument of first exponential function in Eq. (10)
numerator = (x+x.^2)./(1-c*x).^2*(N_i*N3/(Omega_i*Omega_ri*P_ri*P3)) + (x.^2)./(1-c*x).^3 * (b_i^2*Pi/(Omega_i*Omega_ri*P_ri));  %Numerator of the expression in the square roots in Eq. (10)
denominator = 1 + c*x./(1-c*x) * (Pi*Omega_i/(P_ri*Omega_ri)); %Denominator of the expression in the square roots in Eq. (10)

%Compute outage probability using Eq. (10)
outageProbability = 1 - 2* exp(expArgument).*sqrt(numerator./denominator) .* besselk(1,2*sqrt(numerator.*denominator));
if c>0
    outageProbability(x>1/c) = 1;
end

%Compute the integrand of Eq. (16)
integrand = factorIntegrand.*outageProbability;