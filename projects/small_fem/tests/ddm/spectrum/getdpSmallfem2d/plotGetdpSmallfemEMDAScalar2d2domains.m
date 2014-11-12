clear all;
close all;

lAsf = dlmread('lEmdaSf2.csv');
lAgp = dlmread('lEmdaGp2.csv');

figure;
hold on;
plot(real(lAsf), imag(lAsf), 'k+');
plot(real(lAgp), imag(lAgp), 'rx');
hold off;

title('EMDA Scalar 2D -- Getdp vs SmallFem -- Order 1 -- 2 Subdomains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
print('scalEMDAGetdpSmallfem2d2domains.png');
