#!/usr/bin/octave -q

%% Clear
clear all;
close all;

%% Get arguments
args = argv();

%% Load file
load(args{1});

nMesh   = size(data, 1);
nVolume = size(data, 2);
nBorder = size(data, 3);

%% Get Size
it = 0;
for v = 1:nVolume
  for b = 1:nBorder
    if data(1, v, b) ~= 0
      it = it + 1;
    end
  end
end

%% Populate (data)
l2 = zeros(nMesh, it + 1);

for m = 1:nMesh
  it = 1;
  for v = 1:nVolume
    for b = 1:nBorder
      if data(m, v, b) ~= 0
        l2(m, it + 1) = data(m, v, b);
        it = it + 1;
      end
    end
  end
end

%% Populate (mesh size)
for m = 1:nMesh
  l2(m, 1) = 2^(m-1);
end

%% Write
dlmwrite(args{2}, l2, 'delimiter', ',');
