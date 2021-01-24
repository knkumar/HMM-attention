% recalculate summary for HMM
function [time_in_states] = recalcSummary(time_in_states, recalc, foldername, savename)    
    if recalc
        disp(foldername.name);
        load([foldername.folder,'/',foldername.name]);
        [hmm, time_mat] = fix_center_states(hmm,targets);
        time_in_states.(savename) = time_mat;
        [summary,hmm] = generateSummary(hmm, targets);
        save([foldername.folder,'/fixed_',foldername.name],'targets','hand','sim','hmm','summary');
    end
end

