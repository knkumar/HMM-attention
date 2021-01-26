% TODO: check this file for consitency
fnames = dir('../fits/fix*.mat');
data = table();

for id = 1:length(fnames)
	foldername = fnames(id);
	load([foldername.folder, '/', foldername.name]);
	disp(foldername.name)
    
    blocks = summary.respondTargets(1).agg.blocks;
    for block_id = 1:length(blocks)
        curr_block = blocks{block_id};

        data_T = fillData(summary.respondTargets(1).(curr_block), 0, foldername, block_id);
        data_F = fillData(summary.respondFoils(1).(curr_block), 1, foldername, block_id);
        
        if isempty(data)
            data = [data_T; data_F];
        else
            data = [data; data_T; data_F];
        end
   
    end % end block
    writetable(data, [foldername.name(8:8+13), 'data_model.csv']);
    data = table();
end % end subject

% writetable(data, 'data_model_all.csv');


function subject_data = fillData(block_data, foil, foldername, block_id)
    TF = ['T','F'];
    subject_data = table();
    trials = block_data.trialId';
    n = length(trials);
    subject_data.subject = repmat({foldername.name(8:8+13)}, n, 1);
    subject_data.TorF = repmat(TF(foil+1), n, 1);
    
    subject_data.block_id = repmat(block_id, n, 1);
    subject_data.glitches = cellfun(@(x) length(x)-1, block_data.allRes');
    subject_data.first_glitch = cellfun(@(x) x(1), block_data.allRes');
    subject_data.trial_ID = cellfun(@(x) x(3+foil:end), trials, 'UniformOutput', false);
    subject_data.position = cellfun(@(x) x(2+foil), trials);
    subject_data.target_stim = cellfun(@(x) x(1), trials);
    subject_data.stim = cellfun(@(x) x(1+foil), trials);
    subject_data.pres_time = block_data.presTimes';
    subject_data.reached = block_data.reached';
    subject_data.RT = block_data.firstResC';
    subject_data.HT = cellfun(@(x) x(end), block_data.halt_time);
    subject_data.start_state = block_data.start_state';
    subject_data.start_hmm_state = block_data.start_hmm_state';
    % halt, current and stop are hmm states
    subject_data.halt_state = block_data.halt_state';
    subject_data.curr_state = block_data.curr_state';
    subject_data.stop_state = block_data.stop_state';
    subject_data.obj_state = block_data.obj_state';
    subject_data.allRes = block_data.allRes';
end
