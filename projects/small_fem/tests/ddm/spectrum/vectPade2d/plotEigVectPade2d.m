clear all;
close all;

l1 = dlmread('l4P1.csv');
l2 = dlmread('l4P2.csv');
l4 = dlmread('l4P4.csv');
l8 = dlmread('l4P8.csv');

figure;
hold on;
plot(real(l1), imag(l1), 'ko');
plot(real(l2), imag(l2), 'b+');
plot(real(l4), imag(l4), 'rs');
plot(real(l8), imag(l8), 'gx');
hold off;

title('OSRC Vectorial 2D -- Order 4');
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'Pade 1', 'Pade 2', 'Pade 4', 'Pade 8'});

axis square;
print('vectPade2d.png');
