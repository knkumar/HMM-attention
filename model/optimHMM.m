function [cost] = optimHMM(x, targets, traj, sim)
  
  % x vector is organized as follows
  % 4 parameters for velocity
  % 4 parameters for distance - possible to increase to 14
  % 4 for orientation
  % 4 for acceleration
  % 49 parameters for transition
  
  hmm = getHMM_null(x);

  trajopt.data.x = [traj.data.x(1:10000), traj.data.x(end-10000:end) ];
  trajopt.data.y = [traj.data.y(1:10000), traj.data.y(end-10000:end) ];
  trajopt.data.period = [traj.data.period(1:10000), traj.data.period(end-10000:end) ];
  trajopt.data.block = [traj.data.block(1:10000), traj.data.block(end-10000:end) ];
  trajopt.data.trial = [traj.data.trial(1:10000), traj.data.trial(end-10000:end) ];
  trajopt.data.trialID = [traj.data.trialID(1:10000), traj.data.trialID(end-10000:end) ];
  centerState = 0;

  [targets,trajopt,sim,hmm] = analyze10(targets,trajopt,sim,hmm);
  
  optTrial = strcat(trajopt.data.block, trajopt.data.trial);
  
  uniqOptTrial = unique(optTrial);
  cost = 0;
  for i = 1:length(uniqOptTrial)
    trialID = strcmp(optTrial, uniqOptTrial(i));
    trialName = trajopt.data.trial{i};
    foilFlag = isletter(trialName(2));
    foilLoc = -1; targetLoc = -1;
    % if current presentation is a foil use foil location
    if foilFlag
      foilLoc = str2double(trialName(3));
    else
      targetLoc = str2double(trialName(2));
    end

    % circular shift forward by 200ms
    % circshift shifts elements right by index
    if max(find(trialID)) < (length(trialID)-20)
      trialID = circshift(trialID, 20);
    end

    % trialID is shifted forward by 200ms here onward
    startTrialIndex = find(trialID, 1);

    estState = hmm.est.objState(trialID);
    velCurr = hmm.calc.velocity(1,trialID);
    totVel = sum(abs(velCurr));
    % nk 4/7/2019 - calculate first transition 
    firstTrans = find(diff(estState));
    
    % all movement towards target are correct
    correct = (estState == targetLoc);
    firstCorrect = find(correct, 1);
    
    % any movement not towards target or center is incorrect
    incorrect = ~((estState == targetLoc) ); %| (estState == centerState) );
    incorrect = incorrect | (estState == foilLoc);

    % NK: 1/15/2015 - added distance for boostFar                                                 
    dist = targets.calc.distance(trialID,:); 
    
    % find number of transitions to targets and nontargets
    target_trans = length(find(diff(correct)));
    nontarget_trans = length(find(diff(estState(incorrect))));
    % nontarget_states = sum(incorrect(10:end));

    if ~isempty(firstCorrect)
      boostFar = dist(firstCorrect, targetLoc); %encourage farther distances in the model
    else
      boostFar = 900;
    end
    % clip boosting to 900 units distance
    if boostFar > 900
      boostFar = 900;
    end
    boostFar = round(boostFar);
    
    
    mappingBoost = linspace(-2,2,900);
    lagPenalty = startTrialIndex - firstCorrect; % penalize later shifts

    cost_curr = 1 - target_trans - (2*nontarget_trans) + mappingBoost(boostFar);
    cost = cost - cost_curr;
    
    % if ((target_trans+nontarget_trans) == 0) && (totVel > 60) && ~foilFlag
    %    cost = cost + nontarget_states;
    % elseif foilFlag
    %    cost = cost + 0;
    % else
    %    cost = cost - cost_curr;
    % end


  end
  
end
