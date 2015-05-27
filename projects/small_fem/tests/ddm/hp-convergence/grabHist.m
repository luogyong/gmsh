function grabHist(type, ddm, nMsh, nVol, nBrd, file)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Usage: grabJist(type, ddm, nMsh, nVol, nBrd, file)         %%
  %%   type      is 'vector' or 'scalar'.                       %%
  %%   ddm       is the DDM formulation (e.g. osrc, jfl).       %%
  %%   nMsh      is the total number of mehes.                  %%
  %%   nVol      is the total number of volume orders.          %%
  %%   nBrd      is the total number of border orders.          %%
  %%   file      is the name of the output file.                %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% Polpulate data: iterate on Mesh, Volume orders and Border orders
  for m = 1:nMsh
    for v = 1:nVol
      for b = 1:nBrd

        % Insert file
        try
          name          = sprintf('%s_%s_%d_%d_%d.hist', type, ddm, m, v, b);
          hist{m, v, b} = dlmread(name);
        catch
          hist{m, v, b} = 0;
        end
      end
    end
  end

  save(file, 'hist');
end
