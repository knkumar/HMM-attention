function [summary, hmm] = generateSummary(hmm, targets)
  % get summary statistics on the model
  % FCR - reached
  % FCH - halted
  pBlocks = strcat(targets.data.period,targets.data.block);
  [~,idx] = unique(pBlocks);
  uniqBId = pBlocks(sort(idx));
  % strucutre to store responses to targets and foils
  respondFoils = struct(); respondTargets = struct();  
  respondFoils.Red = struct(); respondFoils.Size = struct(); respondFoils.Square = struct();
  respondTargets.Red = struct(); respondTargets.Size = struct(); respondTargets.Square = struct();
  
  statusMat = {'TT','TF';'FT','FF'};
  
  % for each period
  prevTFInd = 0;
  prevResponse = 0;
  reachedTmat = [];
  reachedFmat = [];

  orientedTmat = [];
  orientedFmat = [];
  % for each block calculate summary data
  for p = 1:length(uniqBId)
    periodBlk = uniqBId{p};
    blkId = strcmp(pBlocks,periodBlk);
    blkTrials = targets.data.trial(blkId);
    [~,idx] = unique(targets.data.trial(blkId));
    uniqueTrials = blkTrials(sort(idx));
    lastTrial = 1;
    respondFoils = initResponses(respondFoils, periodBlk);
    respondTargets = initResponses(respondTargets, periodBlk);
    
    blkStates = hmm.est.state(blkId);
    blkObjStates = hmm.est.objState(blkId);
    blkTargetStates = targets.data.state.current(blkId);
    stimDist = targets.calc.distance(blkId,:);
    stimAngles = targets.calc.angles(blkId,:);
    stimVel         = hmm.calc.velocity(blkId);
    blkReachedT = [];
    blkReachedF = [];
    blkOrientedT = [];
    blkOrientedF = [];
    
    % for each trial in block
    for j = 1:length(uniqueTrials)
      currTrial = uniqueTrials{j};
      trialId = strcmp(blkTrials,currTrial);
      
      trialId = circshift(trialId,20);
      %Get States
      states          = blkStates(trialId);
      objStates       = blkObjStates(trialId);
      targetStates    = blkTargetStates(trialId);
      if isletter(currTrial(2))         
        currTrialPos  = str2double(currTrial(3));
      else
        currTrialPos  = str2double(currTrial(2));
      end
      % get angle and distance to stimuli for the current presentation
      currDist        = stimDist(trialId, currTrialPos);
      currAngle       = stimAngles(trialId, currTrialPos);
      currVel         = stimVel(trialId);

      % Generate a vector for last Trial states
      lastTrialMat = ones(1,length(states)) * lastTrial(end);
      % compute changes in states
      
      [changesIdx, stopping] = getTrans(states);
      [OchangesIdx, ignore] = getTrans(objStates);
      
      % if there are any changes before 4 time steps
      if sum(changesIdx < 4)
        spurious.(periodBlk).(currTrial).changes = OchangesIdx;
        spurious.(periodBlk).(currTrial).states = states;
        spurious.(periodBlk).(currTrial).objStates = objStates;
      end
      matches = (objStates==targetStates);
      matchesPrev  = (objStates == lastTrialMat);
      %firstMatch = find(matches);
      % should match the correct state and should have changed
      % state not remained stationary
      firstCorrect = intersect(find(matches),OchangesIdx) + 20;
      
      
      %minDistToTarget = min(currDist);      
      %minAngleToTarget = min(currAngle(currAngle>0));           
      stimStatus = statusMat{prevTFInd+1, isletter(currTrial(2))+1 };
      reached = sum(currDist' < 200 & currVel <10);    
      halting = (currDist' > 50 & currVel<2);
      % time when halted
      timeReached = find(reached,1);      
      notReached = ones(size(halting));
      if ~isempty(timeReached)
        notReached(timeReached:end) = 0;
      end
      halting = halting & notReached; %& stopping;
      halt_time = find(halting);

      if isempty(halt_time)
        halt_time = -10;
      end

      haltId = blkId; haltId(haltId) = trialId;
      haltId(haltId) = halting;
      
      % if sum(diff(halting) == 0) > 4
      %   hmm.est.state(haltId) = 7;
      %   hmm.est.objState(haltId) = 0;
      % end

      if isletter(currTrial(2))
        respondFoils = fillResponses(respondFoils, changesIdx, stimStatus, ...
                            OchangesIdx, matches, matchesPrev, length(states), firstCorrect, ...
                            periodBlk, currTrial, reached, prevResponse, halting, currTrial(2), halt_time);
      else
        respondTargets = fillResponses(respondTargets, changesIdx, stimStatus, ...
                            OchangesIdx, matches, matchesPrev, length(states), firstCorrect, ...
                            periodBlk, currTrial, reached, prevResponse, halting, 0, halt_time);
      end
      
      lastTrial = targetStates(end);
      prevTFInd = isletter(currTrial(2));
      prevResponse = ~(isempty(firstCorrect));
    end
    
    respondTargets = calcAgg(respondTargets, periodBlk, p);
    respondFoils = calcAgg(respondFoils, periodBlk, p);
%     reachedTmat = [reachedTmat; blkReachedT];
%     reachedFmat = [reachedFmat; blkReachedF];
%     orientedFmat = [orientedFmat; blkOrientedF];
%     orientedTmat = [orientedTmat; blkOrientedT];
  end
  
  summary = struct('respondFoils',respondFoils,'respondTargets',respondTargets,'spurious',spurious);
    %'reachedFmat',reachedFmat, 'reachedTmat', reachedTmat, 'orientedTmat', orientedTmat, 'orientedFmat', orientedFmat);
  
end


% function [data] = addResponses(data, reached, changesIdx, stimStatus, OchangesIdx,...
%                     matches,matchesPrev, states, firstCorrect, periodBlk, ...
%                     currTrial, prevResponse, halting)
    
%     data = fillResponses(data, changesIdx, stimStatus, ...
%                 OchangesIdx, matches, matchesPrev, length(states), firstCorrect, ...
%                 periodBlk, currTrial, reached, prevResponse, halting);
% end


function data = calcAgg(resStruct, periodBlk, id)
  data = resStruct;
  for structlen = 1:4
    data(structlen).agg.blocks{id} = periodBlk;
    targetType = periodBlk(4:end-1);
    curr = data(structlen).(periodBlk);
    FR = curr.firstRes(curr.firstRes>0);
    RES = curr.Res(curr.Res<95);    
    NC = curr.Nchanges;
    FC = curr.firstResC(curr.firstResC>0);
    propC = curr.firstResC > 0;
    
    FC = checkEmpty(FC); NC = checkEmpty(NC);
    propC = checkEmpty(propC); RES = checkEmpty(RES);
    FR = checkEmpty(FR);
    
    data(structlen) = fillAggwithID(data(structlen), curr, id, FR, RES, FC, NC, propC, targetType);
    % add this after calcAgg calculations
    data(structlen).aggC = fillConditionals(data(structlen), curr, id);
  end

end

function data_aggC = fillConditionals(data, curr, id)
    % Response conditional on which stimulus was presented before
    % - conditions {FF,FT,TF,TT}
    %
    % fc - first correct
    % fcr - first correct response
    % fch - first correct halted?
    conditions = {'FF','FT','TF','TT'};
    for idx = 1:length(conditions)
        condData = curr.(conditions{idx});

        FC = condData.firstResC; FC = checkEmpty(FC);
        FCR = condData.firstResC((condData.firstResC>0) & (condData.reached>0)); 
        FCR = FCR(FCR > 0); FCR = checkEmpty(FCR);
        FCH = condData.firstResC((condData.firstResC>0) & (condData.halted>0));
        FCH = FCH(FCH > 0); FCH = checkEmpty(FCH);

        propC = FC>0; propC = checkEmpty(propC);
        propCR = (condData.firstResC>0) & (condData.reached>0);
        propCR = checkEmpty(propCR);
        propCH = (condData.firstResC>0) & (condData.halted>0);
        propCH = checkEmpty(propCH);

        data.aggC.(conditions{idx}).meanpropC(id) = mean(propC);        
        data.aggC.(conditions{idx}).meanpropCR(id) = mean(propCR);
        data.aggC.(conditions{idx}).meanpropCH(id) = mean(propCH);
        
        data.aggC.(conditions{idx}).meanFC(id) = mean(FC(FC>0));
        data.aggC.(conditions{idx}).meanFCR(id) = mean(FCR(FCR>0));
        data.aggC.(conditions{idx}).meanFCH(id) = mean(FCH(FCH>0));
        
        FR = condData.firstRes(condData.firstRes>0);
        FR = checkEmpty(FR);
        data.aggC.(conditions{idx}).meanFR(id) = mean(FR(FR>0));
        
        NC = curr.(conditions{idx}).Nchanges;
        NC = checkEmpty(NC);
        data.aggC.(conditions{idx}).meanNR(id) = mean(NC);
        data.aggC.(conditions{idx}).meanNC(id) = mean(propC);
        data.aggC.(conditions{idx}).meanprop(id) = mean(NC>0);
        % number of data points
        data.aggC.(conditions{idx}).num_prop = sum(NC>0);
        data.aggC.(conditions{idx}).num_propC = length(propC);
        data.aggC.(conditions{idx}).num_FC = sum(FC>0);

        data.aggC.(conditions{idx}).errFR(id) = sqrt(var(FR)/length(FR));
        data.aggC.(conditions{idx}).errFC(id) = sqrt(var(FC)/length(FC));
        data.aggC.(conditions{idx}).errFCR(id) = sqrt(var(FCR)/length(FC));
        data.aggC.(conditions{idx}).errFCH(id) = sqrt(var(FCH)/length(FC));
        data.aggC.(conditions{idx}).errNR(id) = sqrt(var(NC)/length(NC));
        data.aggC.(conditions{idx}).errpropC(id) = sqrt(var(propC)/length(propC));
        data.aggC.(conditions{idx}).errpropCR(id) = sqrt(var(propCR)/length(propCR));
        data.aggC.(conditions{idx}).errpropCH(id) = sqrt(var(propCH)/length(propCH));
        %data.aggC.(conditions{idx}).errNC(id) = sqrt(var(propC)/length(propC));
  
    end
    data_aggC = data.aggC;
end

function data = checkEmpty(data)
  if isempty(data)
    data = 0;
  end
end

function data = fillAggwithID(data, curr, id, FR, RES, FC, NC, propC, targetType )
  % fillAggwithID calculates aggregate for each block
  % id = block, FR = first transition, FC = first correct transition
  % NC = Number of correct transitions
  % propC = proportion of correct transition (hit rate)
  

  FCR = curr.firstResC((curr.firstResC>0) & (curr.reached>0));
  FCR = FCR(FCR>0); FCR = checkEmpty(FCR);
  FCH = curr.firstResC((curr.firstResC>0) & (curr.halted>0));
  FCH = FCH(FCH>0); FCH = checkEmpty(FCH);
  propCR = (curr.firstResC>0) & (curr.reached>0);
  propCR = checkEmpty(propCR);
  propCH = (curr.firstResC>0) & (curr.halted>0);
  propCH = checkEmpty(propCH); 
  
  agg = data.agg;
  agg.meanFR(id) = mean(FR); % first response
  agg.meanRes(id) = mean(RES); % all responses
  agg.meanFC(id) = mean(FC(FC>0)); % first correct
  agg.meanFCR(id) = mean(FCR);
  agg.meanFCH(id) = mean(FCH);
  agg.meanNR(id) = mean(NC); % number of transitions

  %agg.meanNC(id) = mean(propC);
  
  agg.meanpropC(id) = mean(propC);
  agg.meanpropCR(id) = mean(propCR);
  agg.meanpropCH(id) = mean(propCH);
  agg.meanprop(id) = mean(NC>0);
  % number of data points
  agg.num_prop(id) = sum(NC>0);
  agg.num_propC(id) = length(propC);
  agg.num_FC(id) = sum(FC>0);
  %disp(data.agg.prop(id))

  % agg.varFR(id) = var(FR);
  % agg.varRes(id) = var(RES); % all responses
  % dagg.varFC(id) = var(FC); % first correct
  % agg.varNR(id) = var(NC); % number of transitions
  
  agg.errFR(id) = sqrt(var(FR)/length(FR));
  agg.errRes(id) = sqrt(var(RES)/length(RES));
  agg.errFC(id) = sqrt(var(FC)/length(FC));
  agg.errFCR(id) = sqrt(var(FCR)/length(FCR));
  agg.errFCH(id) = sqrt(var(FCH)/length(FCH));
  agg.errNR(id) = sqrt(var(NC)/length(NC));
  agg.errpropC(id) = sqrt(var(propC)/length(propC));
  agg.errpropCR(id) = sqrt(var(propCR)/length(propCR));
  agg.errpropCH(id) = sqrt(var(propCH)/length(propCH));
  
  % fill mean for each target type
  if isstruct(data.(targetType)) & isempty(fieldnames(data.(targetType)))
    data.(targetType) = structfun(@(x) x(id), agg, 'UniformOutput',false);
  else
    data.(targetType) = [data.(targetType), structfun(@(x) x(id), agg, 'UniformOutput',false) ];
  end

  data.agg = agg;

end

function responseStruct = initResponses(responseStruct, periodBlk)
  for structlen = 1:4
    responseStruct(structlen).(periodBlk).trialId = {};
    responseStruct(structlen).(periodBlk).Res = [];
    responseStruct(structlen).(periodBlk).firstRes = [];
    responseStruct(structlen).(periodBlk).presTimes = [];
    responseStruct(structlen).(periodBlk).firstResC = [];
    responseStruct(structlen).(periodBlk).prev = [];
    responseStruct(structlen).(periodBlk).Nchanges = [];
    responseStruct(structlen).(periodBlk).Ochanges = [];
    responseStruct(structlen).(periodBlk).reached = [];
    responseStruct(structlen).(periodBlk).halted = [];
    responseStruct(structlen).(periodBlk).halt_time = {};
    responseStruct(structlen).(periodBlk).prevResponse = [];
    responseStruct(structlen).(periodBlk).stim = [];
    
    names = {'FF','FT','TF','TT'};
    for i = 1:length(names)
      responseStruct(structlen).(periodBlk).(names{i}).Ochanges = [];
      responseStruct(structlen).(periodBlk).(names{i}).Nchanges = [];
      responseStruct(structlen).(periodBlk).(names{i}).firstResC = [];
      responseStruct(structlen).(periodBlk).(names{i}).firstRes = [];
      responseStruct(structlen).(periodBlk).(names{i}).Res = [];
      responseStruct(structlen).(periodBlk).(names{i}).reached = [];
      responseStruct(structlen).(periodBlk).(names{i}).halted = [];
      responseStruct(structlen).(periodBlk).(names{i}).prevResponse = [];
    end
  end
end

function timeGrp = getTimeGrp(presTime)
   if presTime < 60
    timeGrp = 1;
  else
    if presTime < 80
      timeGrp = 2;
    else
      timeGrp = 3;
    end
  end
end


function [data] = fillResponses(data, changesIdx, ptInd, OchangesIdx,...
      matches, matchesPrev, presTime, firstCorrect, periodBlk, trial, reached, prevResponse, halting, foilInd, halt_time)
    
  % prev measures the states when moving towards previous state
  
  if isempty(firstCorrect)
    % mark responses for computing aggregates
    firstCorrect = -10;
  end
  
  changeFlag = ones(1,length(matches));
  if OchangesIdx ~= 95
    changeFlag(OchangesIdx) = 1;
  end
  % find changes that are not previous object
  changeFlag = find(changeFlag & (~matchesPrev));
  changes = intersect(changeFlag, OchangesIdx);
  %changes = changesIdx(changesIdx > 4);
  
  if isempty(changes)
    changes = -10;
  end
    
  % 3 groups (<60,<80,>80)
  timeGrp = getTimeGrp(presTime);
  %% set data for each trial
  data(1).(periodBlk).(trial).changes = changesIdx;
  data(1).(periodBlk).(trial).matches = matches;
  data(1).(periodBlk).(trial).matchesprev = matchesPrev;
  data(1).(periodBlk).trialId = [data(1).(periodBlk).trialId, trial];

  data(timeGrp+1).(periodBlk).(trial).changes = changesIdx;
  data(timeGrp+1).(periodBlk).(trial).matches = matches;
  data(timeGrp+1).(periodBlk).(trial).matchesprev = matchesPrev;
  
  % compute statistics for each pBlocks
  function [dat] = computeStatsforDat(dat)
    dat.firstRes = [dat.firstRes, changes(1)]; % first response
    dat.Res = [dat.Res, changesIdx]; % all responses
    dat.presTimes = [dat.presTimes, presTime]; % presentation time
    dat.firstResC = [dat.firstResC, firstCorrect(1)]; % first correct
    dat.reached = [dat.reached, reached];
    dat.halted = [dat.halted, sum(halting)>0];
    dat.halt_time = [dat.halt_time; halt_time];
    if foilInd ~= 0
      dat.stim = [dat.stim,foilInd];
    end
    dat.prevResponse = [dat.prevResponse, prevResponse];
    dat.(ptInd).firstResC = [dat.(ptInd).firstResC, firstCorrect(1)];
    dat.(ptInd).firstRes = [dat.(ptInd).firstRes, changes(1)];
    dat.(ptInd).Res = [dat.(ptInd).Res, changesIdx];    
    dat.(ptInd).reached = [dat.(ptInd).reached, reached];
    dat.(ptInd).halted = [dat.(ptInd).halted, sum(halting)>0];
    dat.(ptInd).prevResponse = [dat.(ptInd).prevResponse,prevResponse];
    if (changes == -10)
        dat.Nchanges = [dat.Nchanges, 0];
        dat.(ptInd).Nchanges = [dat.(ptInd).Nchanges, 0];
    else
        if ~isempty(firstCorrect)
            dat.Nchanges = [dat.Nchanges, length(firstCorrect)]; % number of transitions
            dat.(ptInd).Nchanges = [dat.(ptInd).Nchanges, length(firstCorrect)];
        end
        dat.Ochanges = [dat.Ochanges, changes(1)]; % all responses not towards prev
        dat.(ptInd).Ochanges = [dat.(ptInd).Ochanges, changes(1)];
    end
  end 

  % compute stats for all times
  data(1).(periodBlk) = computeStatsforDat(data(1).(periodBlk) );
  % compute stats conditonal on presentation time
  data(timeGrp+1).(periodBlk) = computeStatsforDat( data(timeGrp+1).(periodBlk) );

end

function [data, stopping] = getTrans(in_data)
  % find transitions in data for each trial
  dData = logical([0 diff(in_data)]);
  stopping = (in_data==2 | in_data==4 | in_data==6);
  data = find(dData);
  % if there are no transiton set it to a large value
  if isempty(data)
    data = 95;
  end
end
