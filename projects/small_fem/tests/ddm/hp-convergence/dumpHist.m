#!/usr/bin/octave -q

%% Clear
clear all;
close all;

%% Get arguments
args = argv();

%% Load file
load(args{1});

nMesh   = size(hist, 1);
nVolume = size(hist, 2);
nBorder = size(hist, 3);

%% Get Size
it = 0;
for v = 1:nVolume
  for b = 1:nBorder
    if length(hist{1, v, b}) ~= 1
      it = it + 1;
    end
  end
end

%% Populate (hist)
dump = zeros(nMesh, it + 1);

for m = 1:nMesh
  it = 1;
  for v = 1:nVolume
    for b = 1:nBorder
      if length(hist{m, v, b}) ~= 1
        dump(m, it + 1) = length(hist{m, v, b});
        it = it + 1;
      end
    end
  end
end

%% Populate (mesh size)
for m = 1:nMesh
  dump(m, 1) = 2^(m-1);
end

%% Write
dlmwrite(args{2}, dump, 'delimiter', ',');
