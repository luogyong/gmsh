function plotIt(test, ddm, type, dom, part, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Usage: plotIt(test, dmm, type, dom, part, k).              %%
%%   test      is the test case (e.g. cylinder, waveguide).   %%
%%   ddm       is the DDM formulation (e.g. osrc4, jfl).      %%
%%   type      is 'vector' or 'scalar'.                       %%
%%   dom       is the number of subdomains.                   %%
%%   part      is the subdomain for L2 error (starting at 1). %%
%%   k         is the wavenumer.                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Capitalize, uppercase and remove numbers
    ctest = capitalize(test);
    ctype = capitalize(type);
    ltype = tolower(type);
    uddm  = toupper(ddm);
    tddm  = regexprep(tolower(ddm), '[0-9]', '');

    %% Root dir
    root = sprintf('%s3DK%dConvergence%s%sTol9Dom%dPart%d', ...
                   ctest, k, uddm, ctype, dom, part);


    %% Read files
    file  = [root, '/', ltype, '_', tddm, '_'];
    f1    = dlmread([file, '1_1.hist']);
    f2    = dlmread([file, '2_2.hist']);
    f3    = dlmread([file, '3_3.hist']);
    f4    = dlmread([file, '4_4.hist']);
    m41   = dlmread([file, '4_1.hist']);
    m42   = dlmread([file, '4_2.hist']);
    m43   = dlmread([file, '4_3.hist']);
    l2Ddm = dlmread([file, '4_4.l2']);
    l2Ddm = l2Ddm(4:end, :);

    %% Name for plots
    name = sprintf('Mixed Order Convergence: %s %s %s -- K: %d', ...
                   ctype, ctest, uddm, k);

    %% Convergence
    figure;
    hold on;

    semilogy([1:size(f1, 1)], f1, '-xk');
    semilogy([1:size(f2, 1)], f2, '-+r');
    semilogy([1:size(f3, 1)], f3, '-ob');
    semilogy([1:size(f4, 1)], f4, '-sg');

    semilogy([1:size(m41, 1)], m41, '-^y');
    semilogy([1:size(m42, 1)], m42, '-vm');
    semilogy([1:size(m43, 1)], m43, '-dc');

    grid;
    hold off;

    title(name);
    xlabel('Iteration');
    ylabel('Residual');

    legend({'Full Order 1',
            'Full Order 2',
            'Full Order 3',
            'Full Order 4',
            'Mixed 4 and 1',
            'Mixed 4 and 2',
            'Mixed 4 and 3'});

    %% Error
    figure;
    hold on;

    semilogy([1:size(l2Ddm, 2)], l2Ddm(1, :), '-xk');
    semilogy([1:size(l2Ddm, 2)], l2Ddm(2, :), '-+r');
    semilogy([1:size(l2Ddm, 2)], l2Ddm(3, :), '-ob');
    semilogy([1:size(l2Ddm, 2)], l2Ddm(4, :), '-sg');

    grid;
    hold off;

    title(name);
    xlabel('Border Order');
    ylabel('L2 Error');

    legend({'Volume Order 1',
            'Volume Order 2',
            'Volume Order 3',
            'Volume Order 4'}, 'location', 'southwest');

    %% Display Error
    disp(l2Ddm);
endfunction

function cap = capitalize(string)
    cap    = lower(string);
    cap(1) = upper(cap(1));
endfunction
