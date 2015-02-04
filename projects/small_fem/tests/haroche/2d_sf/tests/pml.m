%% Clear
close all;
clear all;

%% Problem: 1 / (mu_r * epsilon_r) * A * x - B * lambda * x = 0
%  Solver: pos_gen_non_hermitian (slepc 3.4)
%  Two kinds of PMLs: constant and 1 / x

%% Variation on PML distance
%  |1|x|10|10|10|2|tri|4|4|-|1e-15|

dist_posgen          = [1:5];
dist_posgen_oneOverF = [18879 19123 18898 19188 19008];
dist_posgen_constant = [17644 27758 20506 16009 20231];

figure;
plot(dist_posgen, [dist_posgen_oneOverF; dist_posgen_constant]);
xlabel('PML distance [wavelength]');
ylabel('Lifetime [ms]');
title('Positive General Non-Hermitian (slepc 3.5) -- PML distance');
legend({'PML: 1/x', 'PML: constant'});

%% Variation on PML size
%  |x|5|10|10|10|2|tri|4|4|-|1e-15|

size_posgen          = [1:5];
size_posgen_oneOverF = [19008 18968 18973 18974 18974];
size_posgen_constant = [20231 20231 20231 20231 20231];

figure;
plot(size_posgen, [size_posgen_oneOverF; size_posgen_constant]);
xlabel('PML size [wavelength]');
ylabel('Lifetime [ms]');
title('Positive General Non-Hermitian (slepc 3.5) -- PML size');
legend({'PML: 1/x', 'PML: constant'});

%% Problem: 1 / mu_r * A * x - epsilon_r * B * lambda * x = 0
%  Solver: gen_non_hermitian (slepc 3.4)
%  Two kinds of PMLs: constant and 1 / x

%% Variation on PML distance
%  |1|x|10|10|10|2|tri|4|4|-|1e-15|

dist_gen          = [1:5];
dist_gen_oneOverF = [19124 19200 19097 19211 18974];
dist_gen_constant = [19123 19201 19098 19211 18974];

figure;
plot(dist_gen, [dist_gen_oneOverF; dist_gen_constant]);
xlabel('PML distance [wavelength]');
ylabel('Lifetime [ms]');
title('General Non-Hermitian (slepc 3.5) -- PML distance');
legend({'PML: 1/x', 'PML: constant'});

%% Variation on PML size
%  |x|5|10|10|10|2|tri|4|4|-|1e-15|

size_gen          = [1:5];
size_gen_oneOverF = [18974 18974 18974 18974 18974];
size_gen_constant = [18974 18974 18974 18974 18974];

figure;
plot(size_gen, [size_gen_oneOverF; size_gen_constant]);
xlabel('PML size [wavelength]');
ylabel('Lifetime [ms]');
title('General Non-Hermitian (slepc 3.5) -- PML size');
legend({'PML: 1/x', 'PML: constant'});

%% Comparison between positive and non-positive

figure;
plot(size_gen, [size_posgen_oneOverF; size_gen_oneOverF]);
xlabel('PML size [wavelength]');
ylabel('Lifetime [ms]');
title('PML (1/x) size (slepc 3.5)');
legend({'Positive general non-hermitian', 'General non-hermitian'});

figure;
plot(size_gen, [size_posgen_constant; size_gen_constant]);
xlabel('PML size [wavelength]');
ylabel('Lifetime [ms]');
title('PML (constant) size (slepc 3.5)');
legend({'Positive general non-hermitian', 'General non-hermitian'});

figure;
plot(dist_gen, [dist_posgen_oneOverF; dist_gen_oneOverF]);
xlabel('PML dist [wavelength]');
ylabel('Lifetime [ms]');
title('PML (1/x) distance (slepc 3.5)');
legend({'Positive general non-hermitian', 'General non-hermitian'});

figure;
plot(dist_gen, [dist_posgen_constant; dist_gen_constant]);
xlabel('PML dist [wavelength]');
ylabel('Lifetime [ms]');
title('PML (constant) distance (slepc 3.5)');
legend({'Positive general non-hermitian', 'General non-hermitian'});
