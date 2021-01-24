function mouse_data = normalizeData(data, modRes)
    %%
    % function to sample data to quantify attention in cursor movement
    % 
    %%

    intervals = 0.01:modRes:1.3; % sample every 10 ms
    curr_lower = 0;
    mouse_data = [];
    if size(data,1)  < 1
      return
    end
    
    for i = 1:length(intervals)        
        curr_upper = intervals(i);
        timing = data(:,1);
        timing = timing - timing(1);
        act_data = data(timing >= curr_lower & timing <= curr_upper, :);
        if size(act_data,1) <= 1
          mouse_data = [mouse_data; act_data];
          curr_lower = curr_upper;
          continue;
        end

        meanx = mean(act_data(:,2));
        meany = mean(act_data(:,3));

        % decide on direction of movement
        %diffx = diff(act_data(:,2));
        %diffy = diff(act_data(:,3));
        % Total difference along both x and y directions
        %[valx,idx] = max(abs(diffx));
        %[valy,idy] = max(abs(diffy));
        
        %data_to_store = [curr_upper, act_data(idx,2), act_data(idy,3), ...
        %                 act_data(1,4), act_data(1,5)];
        
        data_to_store = [curr_upper, meanx, meany, act_data(1,4), act_data(1,5)];

        mouse_data = [mouse_data; data_to_store];
        curr_lower = curr_upper;
    end
end