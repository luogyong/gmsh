clear all;
close all;

l0 = dlmread('l0.csv');
l1 = dlmread('l1.csv');
l2 = dlmread('l2.csv');
l3 = dlmread('l3.csv');

figure;
hold on;
plot(real(l0), imag(l0), 'k+');
plot(real(l1), imag(l1), 'ko');
plot(real(l2), imag(l2), 'bo');
plot(real(l3), imag(l3), 'ro');
hold off;

title('OSRC Vectorial 3D');
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'Order 0', 'Order 1', 'Order 2', 'Order 3'});

axis square;
print('vect3d.png');
