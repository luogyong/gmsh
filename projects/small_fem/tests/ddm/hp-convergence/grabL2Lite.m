function grabL2Lite(type, nDom, nMsh, nOrd, file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Usage: grabL2Lite(type, nDom, nMsh, nOrd, file)            %%
%%   type      is 'vector' or 'scalar'.                       %%
%%   nDom      is the total number of subdomains.             %%
%%   nMsh      is the total number of mehes.                  %%
%%   nOrd      is the total number of FEM orders.             %%
%%   file      is the name of the output file.                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

domSum = zeros(1, nDom);
data   = zeros(nMsh, nOrd);

%% Polpulate data: iterate on mesh and FEM orders
for m = 1:nMsh
  for o = 1:nOrd
    % Clear domSum
    domSum = domSum * 0;

    % Iterate on domains
    for d = 0:(nDom - 1)
      try
        name = sprintf('%s_%d_%d_%d.tmp', type, m, o, d);
        domSum(d + 1) = dlmread(name);
      catch
        domSum = domSum * 0;
      end_try_catch
    end

    % Sum domains contribution and store in data
    data(m, o) = sum(domSum);
  end
end

save(file, 'data');
end
