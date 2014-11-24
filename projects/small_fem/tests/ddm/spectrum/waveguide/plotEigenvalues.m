%% Clear
%%%%%%%%
clear all;
close all;

%% Data
%%%%%%%
lsfScalarEmda = dlmread('lsfScalarEmdaWaveguide3dDom2K10O1.csv');
lsfScalarOO2  = dlmread('lsfScalarOO2Waveguide3dDom2K10O1.csv');
lsfScalarOsrc = dlmread('lsfScalarOsrc4Waveguide3dDom2K10O1.csv');

lsfVectorEmda = dlmread('lsfVectorEmdaWaveguide3dDom2K10O0.csv');
lsfVectorJfl  = dlmread('lsfVectorJflWaveguide3dDom2K10O0.csv');
lsfVectorOsrc = dlmread('lsfVectorOsrc4Waveguide3dDom2K10O0.csv');

lgpScalarEmda = dlmread('lgpScalarEmdaWaveguide3dDom2K10O1.csv');
lgpScalarOO2  = dlmread('lgpScalarOO2Waveguide3dDom2K10O1.csv');
lgpScalarOsrc = dlmread('lgpScalarOsrc4Waveguide3dDom2K10O1.csv');

lgpVectorEmda = dlmread('lgpVectorEmdaWaveguide3dDom2K10O0.csv');
lgpVectorJfl  = dlmread('lgpVectorJflWaveguide3dDom2K10O0.csv');
lgpVectorOsrc = dlmread('lgpVectorOsrc4Waveguide3dDom2K10O0.csv');

%% Plot
%%%%%%%

%% Scalar EMDA
figure;
hold on;
plot(real(lsfScalarEmda), imag(lsfScalarEmda), 'k+');
plot(real(lgpScalarEmda), imag(lgpScalarEmda), 'rx');
hold off;
grid;

title('Waveguide 3D: Scalar EMDA -- Getdp vs SmallFem -- Order 0 -- 2 Domains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
axis equal;

%% Scalar OO2
figure;
hold on;
plot(real(lsfScalarOO2), imag(lsfScalarOO2), 'k+');
plot(real(lgpScalarOO2), imag(lgpScalarOO2), 'rx');
hold off;
grid;

title('Waveguide 3D: Scalar OO2 -- Getdp vs SmallFem -- Order 0 -- 2 Domains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
axis equal;

%% Scalar OSRC4
figure;
hold on;
plot(real(lsfScalarOsrc), imag(lsfScalarOsrc), 'k+');
plot(real(lgpScalarOsrc), imag(lgpScalarOsrc), 'rx');
hold off;
grid;

title('Waveguide 3D: Scalar OSRC4 -- Getdp vs SmallFem -- Order 0 -- 2 Domains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
axis equal;

%% Vectorial EMDA
figure;
hold on;
plot(real(lsfVectorEmda), imag(lsfVectorEmda), 'k+');
plot(real(lgpVectorEmda), imag(lgpVectorEmda), 'rx');
hold off;
grid;

title('Waveguide 3D: Vectorial EMDA -- Getdp vs SmallFem -- Order 0 -- 2 Domains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
axis equal;

%% Vectorial JFL
figure;
hold on;
plot(real(lsfVectorJfl), imag(lsfVectorJfl), 'k+');
plot(real(lgpVectorJfl), imag(lgpVectorJfl), 'rx');
hold off;
grid;

title('Waveguide 3D: Vectorial JFL -- Getdp vs SmallFem -- Order 0 -- 2 Domains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
axis equal;

%% Vectorial OSRC
figure;
hold on;
plot(real(lsfVectorOsrc), imag(lsfVectorOsrc), 'k+');
plot(real(lgpVectorOsrc), imag(lgpVectorOsrc), 'rx');
hold off;
grid;

title('Waveguide 3D: Vectorial OSRC4 -- Getdp vs SmallFem -- Order 0 -- 2 Domains -- Wavenumber 10')
xlabel('Real(lambda)');
ylabel('Imag(lambda)');
legend({'SmallFem', 'Getdp'});

xlim([ 0  2]);
ylim([-1 +1]);

axis square;
axis equal;
