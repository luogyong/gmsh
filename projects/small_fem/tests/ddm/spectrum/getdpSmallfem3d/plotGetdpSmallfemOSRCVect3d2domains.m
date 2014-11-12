clear all;
close all;

lsf = dlmread('lOsrcSf2.csv');
lgp = dlmread('lOsrcGp2.csv');

figure;
hold on;
plot(real(lsf), imag(lsf), 'k+');
plot(real(lgp), imag(lgp), 'rx');
hold off;

title('OSRC Vectorial 3D -- Getdp vs SmallFem -- Order 0 -- 2 Subdomains  -- Wavenumber 5')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
print('vectOSRCGetdpSmallfem3d2domains.png');
