clear all;
close all;

l0 = dlmread('lv0osrc.csv');
l1 = dlmread('lv1osrc.csv');
l2 = dlmread('lv2osrc.csv');
l3 = dlmread('lv3osrc.csv');

figure;
hold on;
plot(real(l0), imag(l0), 'm*');
plot(real(l1), imag(l1), 'k+');
plot(real(l2), imag(l2), 'bx');
plot(real(l3), imag(l3), 'ro');
hold off;

title('OSRC Vectorial 3D -- 4 Subdomains -- Wavenumber 5');
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'Order 0', 'Order 1', 'Order 2', 'Order 3'});

axis square;
print('vect3d.png');
