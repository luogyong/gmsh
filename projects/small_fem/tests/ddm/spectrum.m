%% User data %%
root  = 'spectrumVectCir3D';
maxO  = 4;
minO  = 3;
color = ['k' 'r' 'b' 'g' 'c' 'm' 'y'];
mark  = ['+' 'o' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h'];

%% Gather eigenvalues %%
for(i = minO:maxO)
  for(j = minO:i)
    try
      eval(sprintf("l{%d, %d} = dlmread('%s/%d%d.csv');", ...
                   i - minO + 1, j - minO + 1, root, i, j));
      eval(sprintf("name{%d, %d} = 'Volume: %d - Border: %d';", ...
                   i - minO + 1, j - minO + 1, i, j));
    catch
      display(sprintf("Skipping Volume: %d - Border: %d", i, j));
    end_try_catch
  end
end

%% Shifted unit circle %%
x  = [-1:0.01:1];
yU = +sqrt(1 - x.^2);
yD = -sqrt(1 - x.^2);
x  = x + 1;

%% Plot Volume = Border Order %%
figure;
hold on;
for(i = minO:maxO)
  plot(real(l{i - minO + 1, i - minO + 1}),
       imag(l{i - minO + 1, i - minO + 1}), ...
       'linestyle', 'none', ...
       'marker', mark(i - minO + 1), ...
       'color', color(i - minO + 1));
end

for(i = minO:maxO)
  leg{i - minO + 1} = name{i - minO + 1, i - minO + 1};
end
legend(leg);

plot(x, yU, '-k');
plot(x, yD, '-k');
hold off;

%% Plot Volume order = max and Border order = [min, max] %%
figure;
hold on;
for(i = minO:maxO)
  plot(real(l{maxO - minO + 1, i - minO + 1}), ...
       imag(l{maxO - minO + 1, i - minO + 1}), ...
       'linestyle', 'none', ...
       'marker', mark(i - minO + 1), ...
       'color', color(i - minO + 1));
end

for(i = minO:maxO)
  leg{i - minO + 1} = name{maxO - minO + 1, i - minO + 1};
end
legend(leg);

plot(x, yU, '-k');
plot(x, yD, '-k');
hold off;
