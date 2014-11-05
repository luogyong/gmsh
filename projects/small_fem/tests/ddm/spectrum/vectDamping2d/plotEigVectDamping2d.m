clear all;
close all;

l00 = dlmread('l4P4e0.00.csv');
l10 = dlmread('l4P4e0.10.csv');
l20 = dlmread('l4P4e0.20.csv');

figure;
hold on;
plot(real(l00), imag(l00), 'ko');
plot(real(l10), imag(l10), 'b+');
plot(real(l20), imag(l20), 'rs');
hold off;

title('OSRC Vectorial 2D -- Order 4 -- Pade 4');
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'Damping 00%', 'Damping 10%', 'Damping 20%'});

axis square;
print('vectDamping2d.png');
