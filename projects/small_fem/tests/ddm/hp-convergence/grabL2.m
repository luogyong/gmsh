function grabL2(type, ddm, nDom, nMsh, nVol, nBrd, file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Usage: grabL2(type, ddm, nDom, nMsh, nVol, nBrd, file)     %%
  %%   type      is 'vector' or 'scalar'.                       %%
  %%   ddm       is the DDM formulation (e.g. osrc, jfl).       %%
  %%   nDom      is the total number of subdomains.             %%
  %%   nMsh      is the total number of mehes.                  %%
  %%   nVol      is the total number of volume orders.          %%
  %%   nBrd      is the total number of border orders.          %%
  %%   file      is the name of the output file.                %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  domSum = zeros(1, nDom);
  data   = zeros(nMsh, nVol, nBrd);

  %% Polpulate data: iterate on Mesh, Volume orders and Border orders
  for m = 1:nMsh
    for v = 1:nVol
      for b = 1:nBrd
       % Clear domSum
        domSum = domSum * 0;

        % Iterate on domains
        for d = 0:(nDom - 1)
          try
            name = sprintf('%s_%s_%d_%d_%d_%d.tmp', type, ddm, m, v, b, d);
            domSum(d + 1) = dlmread(name);
          catch
            domSum = domSum * 0;
          end_try_catch
        end

        % Sum domains contribution and store in data
        data(m, v, b) = sum(domSum);
      end
    end
  end

  save(file, 'data');
end
