function [hand, targets] = getDataAnalysis(sName, Resolution)
    load(sName);
    % generate mouse data at a specified resolution
    mouse_data = generateMouseData(block_data, 0.01);
    names = fieldnames(mouse_data);
    % generate target and trajectory data
    for i = 1:length(names)
        name = names{i};
        blkTrial = strsplit(name,'_');
        mouse_matrix = mouse_data.(name);
        n = size(mouse_matrix,1);
        
        trial = blkTrial{3};
        if isletter(trial(2))
            current = str2double(trial(3));
            tflag = 0;
        else
            current = str2double(trial(2));
            tflag = 1;
        end
        if isletter(trial(2))
	        tId = str2double(trial(4:end));
        else
	        tId = str2double(trial(3:end));
        end
        
        if i==1
            hand.data.x = mouse_matrix(:,1)';
            hand.data.y = mouse_matrix(:,2)';
            hand.data.period = repmat(blkTrial(1),1,n);
            hand.data.block = repmat(blkTrial(2),1,n);
            hand.data.trial = repmat(blkTrial(3),1,n);
            hand.data.trialID = repmat(tId,1,n);
            %hand.data.pTime = structfun(@(x) x, block_data.);
            targets.data.state.current = repmat(current,1,n);
            targets.data.state.last = NaN(1,n);
            targets.data.tFlag = repmat(tflag,1,n);
        else
            hand.data.x = cat(2,hand.data.x,mouse_matrix(:,1)');
            hand.data.y = cat(2,hand.data.y,mouse_matrix(:,2)');
            hand.data.period = cat(2,hand.data.period, repmat(blkTrial(1),1,n));
            hand.data.block = cat(2,hand.data.block,repmat(blkTrial(2),1,n));
            hand.data.trial = cat(2,hand.data.trial,repmat(blkTrial(3),1,n));
            targets.data.state.current = cat(2, targets.data.state.current, repmat(current,1,n));
            targets.data.state.last = cat(2, targets.data.state.last, repmat(last,1,n));
            targets.data.tFlag = cat(2,targets.data.tFlag, repmat(tflag,1,n));
            hand.data.trialID = cat(2, hand.data.trialID, repmat(tId,1,n));
        end
        last = current;
    end
    targets.tLocs = getPositions(Resolution,1,[50 170 290]);
    targets.Resolution = Resolution;
    targets.width = 160;
    targets.data.period = hand.data.period;
    targets.data.block = hand.data.block;
    targets.data.trial = hand.data.trial;
end
