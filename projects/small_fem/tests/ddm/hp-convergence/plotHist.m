function plotHist(filename)
  %% Static
  color = ['b', 'g', 'r', 'k', 'c', 'm'];
  mark  = ['o', 'x', 's', 'd', 'v'];

  %% Load file
  load(filename);

  nMesh   = size(hist, 1);
  nVolume = size(hist, 2);
  nBorder = size(hist, 3);

  %% Populate mesh vector
  msh = zeros(nMesh, 1);
  for m = 1:nMesh
    msh(m) = 2^(m-1);
  end

  %% Populate iteration vector
  histSize = zeros(nMesh, nVolume, nBorder);

  for m = 1:nMesh
    for v = 1:nVolume
      for b = 1:nBorder
        if length(hist{m, v, b}) == 1
          histSize(m, v, b) = NaN;
        else
          histSize(m, v, b) = length(hist{m, v, b});
        end
      end
    end
  end

  figure;
  hold on;

  it = 1;
  for v = 1:nVolume
    %b = v;
    for b = 1:nBorder
      if isnan(histSize(1, v, b)) == 0
        semilogx(msh, histSize(:, v, b), sprintf('-%s%s', color(v), mark(b)));
        myLegend{it} = sprintf('Volume: %d, Border: %d', v, b);
        it = it + 1;
      end
    end
  end

  title('DDM Convergence')
  xlabel('Mesh size [per wavelength]');
  ylabel('Iteration number [-]');
  legend(myLegend);
  hold off;
  grid;
end
