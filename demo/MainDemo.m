function MainDemo(small)
% script to run all the Demos in sequence
% Input logical variables
%   small - true or false (or 1 or 0) to run with the small example DEM (true)
%       or the full DEM (1

% Change these two lines to align with folder structure in which script is
% run. The resultFolder is needed only if the logical variable printFigs is
% set to true (or 1), for running the demos in an environment without a
% display.
folder = '..\demodata';
resultFolder = '..\results';
printFigs = false; % set to true to print the figures as png in the resultFolder

% read the elevation data in the example
if verLessThan('map','5.0')
    if small
        [Z,R] = geotiffread(fullfile(folder,'SmallDEM.tif')); %#ok<GTRED>
    else
        unzip(fullfile(folder,'n32_e077_1arc_v3.zip'),folder);
        [Z,R] = geotiffread(fullfile(folder,'n32_e077_1arc_v3.tif')); %#ok<GTRED>
    end
else
    if small
        [Z,R] = readgeoraster(fullfile(folder,'SmallDEM.tif'));
    else
        unzip(fullfile(folder,'n32_e077_1arc_v3.zip'),folder);
        [Z,R] = readgeoraster(fullfile(folder,'n32_e077_1arc_v3.tif'));
    end
end

% print the overall message
fprintf('This Demo reproduces Figures 1 to 5 for the paper in IEEE Geoscience and Remote Sensing Letters\n')
if small
    fprintf('(but for a subset of the area)\n')
end

% data for Demo_profile are in a .mat file
m = matfile(fullfile(folder,'profileDemo.mat'));

elapsed = Demo_imageDEM(Z,R,printFigs,resultFolder);
fprintf('Elapsed time = %f seconds\n\n',elapsed);

elapsed = Demo_profile(m,printFigs,resultFolder);
fprintf('Elapsed time = %f seconds\n\n',elapsed);

elapsed = Demo_rotate(Z,printFigs,resultFolder);
fprintf('Elapsed time = %f seconds\n\n',elapsed);

delete(gcp('nocreate'));
poolobj = parpool;
useParallel = poolobj.NumWorkers>1;
elapsed = Demo_shade(Z,R,printFigs,resultFolder,useParallel);
if useParallel
    fprintf('Elapsed time = %f seconds using %d cores\n\n',elapsed,poolobj.NumWorkers);
else
    fprintf('Elapsed time = %f seconds\n\n',elapsed);
end

elapsed = Demo_viewFactor(Z,R,printFigs,resultFolder,useParallel);
if useParallel
    fprintf('Elapsed time = %f seconds using %d cores\n\n',elapsed,poolobj.NumWorkers);
else
    fprintf('Elapsed time = %f seconds\n\n',elapsed);
end
end