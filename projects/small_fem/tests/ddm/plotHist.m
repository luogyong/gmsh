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
    msh(m) = 2^m;
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
  hold on
  for v = 1:nVolume
    %b = v;
    for b = 1:nBorder
      semilogx(msh, histSize(:, v, b), sprintf('-%s%s', color(v), mark(b)));
    end
  end
  hold off
  grid
end
