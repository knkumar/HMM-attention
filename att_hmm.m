function att_hmm(datadir)
    % This is best version of parameters and fits for the HMM
    close all; clc;
    root = pwd; % set root path for model

    addpath(genpath([root,'/applyhatch_pluscolor_bundle']));
    addpath(genpath([root '/pre']));
    addpath(genpath([root '/model']));
    addpath(genpath([root '/post']));

    Resolution = [1200 1600];
    iterations = 200;
    sNames = dir(datadir);
    
    for i = 1:size(sNames,1)
      sName = sNames(i,:);
      if strcmp(sName.name,'.') || strcmp(sName.name,'..')
        continue
      end
      disp(sName.name);
      fitmodel(sName, Resolution, iterations);
    end
    
cd('./fix');
recalc_main    
end


