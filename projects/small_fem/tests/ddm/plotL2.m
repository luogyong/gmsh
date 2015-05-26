function plotL2(filename)
  %% Static
  color = ['b', 'g', 'r', 'k', 'c', 'm'];
  mark  = ['o', 'x', 's', 'd', 'v'];

  %% Load file
  load(filename);

  nMesh   = size(data, 1);
  nVolume = size(data, 2);
  nBorder = size(data, 3);

  %% Populate mesh vector
  msh = zeros(nMesh, 1);
  for m = 1:nMesh
    msh(m) = 2^m;
  end

  %% Plot
  figure;
  hold on;
  for v = 1:nVolume
    for b = 1:nBorder
      loglog(msh, data(:, v, b), sprintf('-%s%s', color(v), mark(b)));
    end
  end
  hold off;
  grid;
end
