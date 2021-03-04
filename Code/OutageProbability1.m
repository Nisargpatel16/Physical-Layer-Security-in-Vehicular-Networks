close all;
gamma_m = (-10:30);
gamma_w_1 = 20;
gamma_w_2 = 10;
gamma_w_3 = 0;
gamma_w_4 = -10;

snr_m = 10.^(gamma_m/10);   
snr_w = 10.^(gamma_w_1/10);

rs = 0.1;

sop_1 = 1 - ((snr_m)./(snr_m + power(2,rs)*snr_w)).*exp((-power(2,rs)+1)./snr_m);
snr_w = 10.^(gamma_w_2/10);
sop_2 = 1 - ((snr_m)./(snr_m + power(2,rs)*snr_w)).*exp((-power(2,rs)+1)./snr_m);
snr_w = 10.^(gamma_w_3/10);
sop_3 = 1 - ((snr_m)./(snr_m + power(2,rs)*snr_w)).*exp((-power(2,rs)+1)./snr_m);
snr_w = 10.^(gamma_w_4/10);
sop_4 = 1 - ((snr_m)./(snr_m + power(2,rs)*snr_w)).*exp((-power(2,rs)+1)./snr_m);


close all;
figure;
semilogy(gamma_m,sop_1,'b*-')
hold on
semilogy(gamma_m,sop_2,'r*-')
hold on
semilogy(gamma_m,sop_3,'g*-')
hold on
semilogy(gamma_m,sop_4,'y*-')
hold off

axis([-10 30 10^-5 1]);
title('Outage probability against received SNR for case 2','FontWeight','Normal');
legend('\Gamma(w)=20 db', '\Gamma=10 db', '\Gamma=0 db', '\Gamma=-10 db');
xlabel('SNR (dB)');
ylabel('Outage Probability');
