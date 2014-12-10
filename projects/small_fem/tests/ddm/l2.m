#!/usr/bin/octave -q

%% Clear
%%%%%%%%
clear all;
close all;

%% Read
%%%%%%%
args = argv();
ddm  = dlmread(args{1});
ana  = dlmread(args{2});

%% Remove header
%%%%%%%%%%%%%%%%
ddm = ddm(4:end, :);
ana = ana(4:end, :);

%% Check
%%%%%%%%
chk = ddm(:, 1) - ana(:, 1);
if(norm(chk) ~= 0)
    disp("Error: matrices don't match");
end

clear chk;

%% Remove point index
%%%%%%%%%%%%%%%%%%%%%
ddm = ddm(:, 2:end);
ana = ana(:, 2:end);

%% Error & norm
%%%%%%%%%%%%%%%
err    = ddm - ana;
rmsErr = norm(norm(err, 'rows')) / norm(norm(ana, 'rows'));

%% Disp
%%%%%%%
format short e;
disp(rmsErr);