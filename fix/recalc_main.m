addpath('../post');
fnames = dir('../fits/S*.mat');
time_in_states = struct();
for id = 1:length(fnames)
	fname = fnames(id);
	model = fname.name;
    savename = strrep(strrep(model,'.mat',''),'-','_');
	sNames = dir([fname.folder,'/',model]);
    time_in_states = recalcSummary(time_in_states,true, fname, savename);
    save('../fits/time_in_states.mat','time_in_states');
end
getOverallStats(model, fname.folder, 'fixed*');