function fitmodel(sName, Resolution, maxiter)
    % get data and experimental event details
    [hand, targets] = getDataAnalysis([sName.folder,'/',sName.name], Resolution);
    sim.screenSize = [1600,1200];
    disp('Analyzing data ...');
    % opitimize parameters using simplex search
    options = optimset('MaxIter', maxiter, 'Display','iter','TolX',1e-1);
    % parameters in model to optimize
    x = [.7, .1, ... % transition
        0, 1, 1, 1, ... % distance
        1, 1, ... % orientation (1, 1, 0, 1,)
        1, 1, 0.2, ... % speed
        1, 1]; % acceleration

    % optimize for the parameters using Nelder Mead simplex
    fval=0;
    [xfin, fval] = fminsearch(@optimHMM,x,options,targets,hand,sim);
    disp(sum(xfin-x))
    hmm = getHMM_null(xfin);
x
    [targets,hand,sim,hmm] = analyze10(targets,hand,sim,hmm);
    targets.data.trial = hand.data.trial;
    % generate summary statistics from fit
    % summary = generateSummary(hmm,targets);
    plot_name = strrep(strrep(sName.name,'.mat',''), '_', '-');
    % name to save file
    % parsave(plot_name,targets,hand,sim,hmm,summary);
    % plot model with fitted parameters
    plotModel(targets, hmm, plot_name, [800, 600], plot_name, ['plots/',plot_name]);
end

function parsave(sName,targets,hand,sim,hmm,summary)
  save(strcat('fits/',sName,'-center.mat'),'targets','hand','sim','hmm','summary');
end
