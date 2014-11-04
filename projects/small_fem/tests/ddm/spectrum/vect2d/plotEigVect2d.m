clear all;
close all;

l0 = dlmread('l0.csv');
l1 = dlmread('l1.csv');
l2 = dlmread('l2.csv');
l3 = dlmread('l3.csv');
l4 = dlmread('l4.csv');
l5 = dlmread('l5.csv');

figure;
hold on;
plot(real(l0), imag(l0), 'k+');
plot(real(l1), imag(l1), 'ko');
plot(real(l2), imag(l2), 'bo');
plot(real(l3), imag(l3), 'ro');
plot(real(l4), imag(l4), 'go');
plot(real(l5), imag(l5), 'yo');
hold off;

title('OSRC Vectorial 2D');
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'Order 0', 'Order 1', 'Order 2', 'Order 3', 'Order 4', 'Order 5'});

axis square;
print('vect2d.png');
