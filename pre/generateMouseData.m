function [mouse_data] = generateMouseData(data, modRes)
    %%
    % function to create a structure for mouse data to measure attention effectively 
    % mouse_data['period''block',trial] = [X,Y]
    %%
    % names of blocks in mouseData
    names = fieldnames(data);
    mouse_data = struct();

    % loop for each session
    for i = 1:length(names)
        blockData = data.(names{i});
        bnames = fieldnames(blockData); % block name
        % loop for each block
        for j = 1:length(bnames)
            trialData = blockData.(bnames{j}).mouse;
            trialNames = fieldnames(trialData);
            % loop for each trial
            for k = 1:length(trialNames)
                trialName = trialNames{k};
                mouse_data_all = trialData.(trialName);
                mouse_data_all(mouse_data_all(:,1) == 0,:) = [];
                if size(mouse_data_all,1) == 0
                  continue;
                end
                % Normalize data for modRes
                mouse_data_block = normalizeData(mouse_data_all, modRes);
                % get the X and Y cursor positions only from data
                mouse_data.([names{i} '_' bnames{j} '_' trialName]) = mouse_data_block(:,2:3);
            end
        end
    end
end