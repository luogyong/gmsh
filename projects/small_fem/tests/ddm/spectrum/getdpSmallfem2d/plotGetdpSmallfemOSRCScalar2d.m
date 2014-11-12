clear all;
close all;

lAsf = dlmread('lOsrcSf.csv');
lAgp = dlmread('lOsrcGp.csv');

figure;
hold on;
plot(real(lAsf), imag(lAsf), 'k+');
plot(real(lAgp), imag(lAgp), 'rx');
hold off;

title('OSRC Scalar 2D -- Getdp vs SmallFem -- Order 1 -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
print('scalOSRCGetdpSmallfem2d.png');
