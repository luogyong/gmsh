function plotL2Lite(filename)
  %% Static
  color = ['b', 'g', 'r', 'k', 'c', 'm'];

  %% Load file
  load(filename);
  nMesh  = size(data, 1);
  nOrder = size(data, 2);

  %% Populate mesh vector
  msh = zeros(nMesh, 1);
  for m = 1:nMesh
    msh(m) = 2^m;
  end

  %% Plot
  figure;
  hold on
  for o = 1:nOrder
    loglog(msh, data(:, o), sprintf('-o%s', color(o)));
  end
  hold off

  %% Slope
  delta = zeros(nOrder, nMesh - 1);

  for m = 1:(nMesh - 1)
    delta(:, m) = ...
    (log10(data(m + 1, :)) - log10(data(m, :))) / ...
    (log10(msh(m + 1))     - log10(msh(m)));
  end

  %% Display
  disp(delta);

end
