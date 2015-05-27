function slope = plotL2(filename)
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
    msh(m) = 2^(m-1);
  end

  %% Plot
  figure;
  hold on;

  it = 1;
  for v = 1:nVolume
    %b = v;
    for b = 1:nBorder
      if data(1, v, b) ~= 0
        loglog(msh, data(:, v, b), sprintf('-%s%s', color(v), mark(b)));
        myLegend{it} = sprintf('Volume: %d, Border: %d', v, b);
        it = it + 1;
      end
    end
  end

  title('L2 Error: Analytic vs FEM')
  xlabel('Mesh size [per wavelength]');
  ylabel('L2 Error [-]');
  legend(myLegend);
  hold off;
  grid;

  %% Grab only full order data
  fullO = zeros(nMesh, nVolume);
  for m = 1:nMesh
    for v = 1:nVolume
      fullO(m, v) = data(m, v, v);
    end
  end

  %% Slopes
  slope = zeros(nVolume, nMesh - 1);

  for m = 1:(nMesh - 1)
    slope(:, m) = (log10(fullO(m + 1, :)) - log10(fullO(m, :))) / ...
                  (log10(msh(m + 1))      - log10(msh(m)));
  end
end
