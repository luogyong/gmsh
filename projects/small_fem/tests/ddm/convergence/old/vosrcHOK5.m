%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% High Order Vectorial OSRC           %%
%% 4 Domains -- Pade 4 -- Wavenumber 5 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

o0 = dlmread('convCylinderOSRC4Dom4K5O0');
o1 = dlmread('convCylinderOSRC4Dom4K5O1');
o2 = dlmread('convCylinderOSRC4Dom4K5O2');
o3 = dlmread('convCylinderOSRC4Dom4K5O3');
o4 = dlmread('convCylinderOSRC4Dom4K5O4');
o5 = dlmread('convCylinderOSRC4Dom4K5O5');

figure;
hold on;
semilogy([1:size(o0, 1)], o0, '-xk', 'linewidth', 3);
semilogy([1:size(o1, 1)], o1, '-+r', 'linewidth', 3);
semilogy([1:size(o2, 1)], o2, '-ob', 'linewidth', 3);
semilogy([1:size(o3, 1)], o3, '-sg', 'linewidth', 3);
semilogy([1:size(o4, 1)], o4, '-dy', 'linewidth', 3);
semilogy([1:size(o5, 1)], o5, '-^m', 'linewidth', 3);
hold off;

legend({'Order 0', 'Order 1', 'Order 2', 'Order 3', 'Order 4', 'Order 5'});
title('High Order Vectorial OSRC (Pade 4): 4 domains -- Wavenumber 5');
xlabel('Iteration');
ylabel('Residual');

print('vosrcHOK5.png');