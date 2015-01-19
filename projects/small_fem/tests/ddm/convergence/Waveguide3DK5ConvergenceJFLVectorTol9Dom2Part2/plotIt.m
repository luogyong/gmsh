clear all;
close all;

f1 = dlmread('vector_jfl_1_1.hist');
f2 = dlmread('vector_jfl_2_2.hist');
f3 = dlmread('vector_jfl_3_3.hist');
f4 = dlmread('vector_jfl_4_4.hist');

m41 = dlmread('vector_jfl_4_1.hist');
m42 = dlmread('vector_jfl_4_2.hist');
m43 = dlmread('vector_jfl_4_3.hist');

%m41 = dlmread('vector_jfl_2_1.hist');
%m42 = dlmread('vector_jfl_3_1.hist');
%m43 = dlmread('vector_jfl_4_1.hist');

l2Ddm = dlmread('vector_jfl_4_4.l2');
l2Ddm = l2Ddm(4:end, :);

%%
figure;
hold on;

semilogy([1:size(f1, 1)], f1, '-xk');
semilogy([1:size(f2, 1)], f2, '-+r');
semilogy([1:size(f3, 1)], f3, '-ob');
semilogy([1:size(f4, 1)], f4, '-sg');

semilogy([1:size(m41, 1)], m41, '-^y');
semilogy([1:size(m42, 1)], m42, '-vm');
semilogy([1:size(m43, 1)], m43, '-dc');

grid;
hold off;

title('Mixed order convergence: vectorial waveguide JFL -- K: 5');
xlabel('Iteration');
ylabel('Residual');

legend({'Full Order 1',
        'Full Order 2',
        'Full Order 3',
        'Full Order 4',
        'Mixed 4 and 1',
        'Mixed 4 and 2',
        'Mixed 4 and 3',
       });

%%
figure;
hold on;

semilogy([1:size(l2Ddm, 2)], l2Ddm(1, :), '-xk');
semilogy([1:size(l2Ddm, 2)], l2Ddm(2, :), '-+r');
semilogy([1:size(l2Ddm, 2)], l2Ddm(3, :), '-ob');
semilogy([1:size(l2Ddm, 2)], l2Ddm(4, :), '-sg');

grid;
hold off;

title('Mixed order convergence: vectorial waveguide JFL -- K: 5');
xlabel('Border Order');
ylabel('L2 Error');

legend({'Volume Order 1',
        'Volume Order 2',
        'Volume Order 3',
        'Volume Order 4'
       }, 'location', 'southwest');