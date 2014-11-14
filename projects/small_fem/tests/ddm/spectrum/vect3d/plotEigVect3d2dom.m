clear all;
close all;

l0 = dlmread('lv0osrc2dom.csv');
l1 = dlmread('lv1osrc2dom.csv');
l2 = dlmread('lv2osrc2dom.csv');
l3 = dlmread('lv3osrc2dom.csv');
l4 = dlmread('lv4osrc2dom.csv');
l5 = dlmread('lv5osrc2dom.csv');

figure;
hold on;
plot(real(l0), imag(l0), 'm*');
plot(real(l1), imag(l1), 'k+');
plot(real(l2), imag(l2), 'bx');
plot(real(l3), imag(l3), 'ro');
plot(real(l4), imag(l4), 'gs');
plot(real(l5), imag(l5), 'yd');
hold off;

title('OSRC Vectorial 3D -- 2 Subdomains -- Wavenumber 5');
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'Order 0', 'Order 1', 'Order 2', 'Order 3', 'Order 4', 'Order 5'});

axis square;
print('vect3d2dom.png');
