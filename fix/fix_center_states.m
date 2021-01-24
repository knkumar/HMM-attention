function [hmm, time_mat] = fix_center_states(hmm, targets)
    % concatenate period and block for unique block identifier
    pBlocks = strcat(targets.data.period,targets.data.block);
    [~,idx] = unique(pBlocks);
    uniqBId = pBlocks(sort(idx));
    time_mat = [];

    % for each period and block calculate states
    for p = 1:length(uniqBId)
        periodBlk = uniqBId{p};
        blkId = strcmp(pBlocks,periodBlk);
        blkTrials = targets.data.trial(blkId);
        [~,idx] = unique(targets.data.trial(blkId));
        uniqueTrials = blkTrials(sort(idx));
        blkStates = hmm.est.state(blkId);

        % for each trial in block
        for j = 1:length(uniqueTrials)
            currTrial = uniqueTrials{j};
            trialId   = strcmp(blkTrials,currTrial);
            states          = blkStates(trialId);
            % remove short durations in a state
            states = remove_small(states, 4);
            all_changes = find(diff(states))+1;
            time_in_state = diff(all_changes);

            % find center states
            center_states = find(states == 7);
            changes_center = intersect(center_states,all_changes);
            if isempty(changes_center) & ~isempty(center_states)
                changes_center = center_states(1);
            end
            time_mat = cat(2,time_mat,time_in_state);

            % correct/handle center states
            for cs = 1:length(changes_center)
                % if next and prev state change is not a center state change center states to next or prev state
                % change to next state take precedence ,i.e.,
                % if next state is non-center change state to next
                % if next state is center change state to prev
                curr_cs = changes_center(cs);
                % find next state and prev state changes
                next = find(all_changes==curr_cs)+1;
                prev = find(all_changes==curr_cs)-1;
                if isempty(next) & curr_cs == 1
                    next = 1;
                end
                if prev==0 | (isempty(prev) & curr_cs == 1)
                    prev = 1;
                end
                if next > length(all_changes)
                    continue;
                end
                next_state = states(all_changes(next));
                prev_state = states(all_changes(prev));

                if next_state < 7
                    states(center_states(center_states>=curr_cs & center_states<all_changes(next) )) = next_state;
                end
                if prev_state < 7
                    states(center_states(center_states>=curr_cs& center_states>all_changes(prev))) = prev_state;
                end
            end
            blkStates(trialId) = states;
        end
        blkStates = fix_post(blkStates,7);
        blkStates = fix_post(blkStates,8);
        blkStates = remove_small(blkStates,5);

        hmm.est.state(blkId) = blkStates;
        hmm.est.objState(hmm.est.state == 1 | hmm.est.state == 2) = 1;
        hmm.est.objState(hmm.est.state == 3 | hmm.est.state == 4) = 2;
        hmm.est.objState(hmm.est.state == 5 | hmm.est.state == 6) = 3;
        hmm.est.objState(hmm.est.state == 7) = 0;
        hmm.est.objState(hmm.est.state == 8) = 4;
    end
end

function [states] = remove_small(states, window_size)
    % function to remove states which are smaller than window_size
    all_changes = find(diff(states))+1;
    time_in_state = diff(all_changes);
    unintended_jump = [time_in_state < window_size];
    %disp(size(unintended_jump))
    changes_to_remove = all_changes(unintended_jump);
    times_to_remove = time_in_state(unintended_jump);
    % remove short durations in a state
    for remove = 1:length(changes_to_remove)
        remove_id = changes_to_remove(remove);
        remove_time = times_to_remove(remove);
        states(remove_id:(remove_id+remove_time-1)) = states(remove_id-1);
    end
end

function states = fix_post(states, num)
    num_states = find(states==num);
    for stateid = 1:length(num_states)
        curr = num_states(stateid);
        if curr > 2
            states(curr) = states(curr-1);
        end
    end
end