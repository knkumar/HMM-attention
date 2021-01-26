function [currentState, vStore, logP] = hmmviterbi_multi2(seq,tr,e,trialId, trial,varargin)
%HMMVITERBI calculates the most probable state path for a sequence.
%   STATES = HMMVITERBI(SEQ,TRANSITIONS,EMISSIONS) given a sequence, SEQ,
%   calculates the most likely path through the Hidden Markov Model
%   specified by transition probability matrix, TRANSITIONS, and emission
%   probability matrix, EMISSIONS. TRANSITIONS(I,J) is the probability of
%   transition from state I to state J. EMISSIONS(K,L) is the probability
%   that symbol L is emitted from state K.
%
%   HMMVITERBI(...,'SYMBOLS',SYMBOLS) allows you to specify the symbols
%   that are emitted. SYMBOLS can be a numeric array or a cell array of the
%   names of the symbols.  The default symbols are integers 1 through N,
%   where N is the number of possible emissions.
%
%   HMMVITERBI(...,'STATENAMES',STATENAMES) allows you to specify the
%   names of the states. STATENAMES can be a numeric array or a cell array
%   of the names of the states. The default statenames are 1 through M,
%   where M is the number of states.
%
%   This function always starts the model in state 1 and then makes a
%   transition to the first step using the probabilities in the first row
%   of the transition matrix. So in the example given below, the first
%   element of the output states will be 1 with probability 0.95 and 2 with
%   probability .05.
%
%   Examples:
%
%     tr = [0.95,0.05;
%             0.10,0.90];
%
%     e = [1/6,  1/6,  1/6,  1/6,  1/6,  1/6;
%            1/10, 1/10, 1/10, 1/10, 1/10, 1/2;];
%
%       [seq, states] = hmmgenerate(100,tr,e);
%       estimatedStates = hmmviterbi(seq,tr,e);
%
%       [seq, states] = hmmgenerate(100,tr,e,'Statenames',{'fair';'loaded'});
%       estimatesStates = hmmviterbi(seq,tr,e,'Statenames',{'fair';'loaded'});
%
%   See also HMMGENERATE, HMMDECODE, HMMESTIMATE, HMMTRAIN.

%   Reference: Biological Sequence Analysis, Durbin, Eddy, Krogh, and
%   Mitchison, Cambridge University Press, 1998.

%   Copyright 1993-2011 The MathWorks, Inc.


% tr must be square

numStates = size(tr,1);
checkTr = size(tr,2);
if checkTr ~= numStates
  error(message('stats:hmmviterbi:BadTransitions'));
end

% number of rows of e must be same as number of states

checkE = size(e,1);
if checkE ~= numStates
  error(message('stats:hmmviterbi:InputSizeMismatch'));
end

numEmissions = size(e,2);
customStatenames = false;

% deal with options
if nargin > 3
  okargs = {'symbols','statenames'};
  [symbols,statenames] = ...
    internal.stats.parseArgs(okargs, {[] []}, varargin{:});
  
  if ~isempty(symbols)
    numSymbolNames = numel(symbols);
    if ~isvector(symbols) || numSymbolNames ~= numEmissions
      error(message('stats:hmmviterbi:BadSymbols'));
    end
    [~, seq]  = ismember(seq,symbols);
    if any(seq(:)==0)
      error(message('stats:hmmviterbi:MissingSymbol'));
    end
  end
  if ~isempty(statenames)
    numStateNames = length(statenames);
    if numStateNames ~= numStates
      error(message('stats:hmmviterbi:BadStateNames'));
    end
    customStatenames = true;
  end
end


% work in log space to avoid numerical issues
L = length(seq);
if any(seq(:)<1) || any(seq(:)~=round(seq(:))) || any(seq(:)>numEmissions)
  error(message('stats:hmmviterbi:BadSequence', numEmissions));
end
currentState = zeros(1,L);
if L == 0
  return
end
tr(tr==0) = 1e-100; %NK : 12/28 processing here 
logTR = log(tr);
e(e==0) = 1e-100;
logE = log(e);

% allocate space
pTR = zeros(numStates,L);
% assumption is that model is in state 1 at step 0
v = -Inf(numStates,1);
v(1,1) = 0; % probability of 1 
vOld = v;

trialChanges = diff(trialId);
vStore = nan(size(pTR));

% loop through the model
% [NK - 6/3/2017 : Create a variable to keep track of the number of times
% states are changed for a given trial
boost_count = 0;
dCount = 1;

% for every observed sequence
for count = 1:L 
  currTrial = trial{count};
  if (count < L) && trialChanges(count)
    dCount = 1;
  else
    dCount = dCount + 1;
  end
  
  if isletter(currTrial(2))
    targetState = [0,0];
  else
    targetPos = str2double(currTrial(2));
    targetState = [(targetPos*2)-1, targetPos*2];
  end
  boost_count = zeros(1,numStates);
  
  % for every state
  for state = 1:numStates 
    % for each state we calculate
    % v(state) = e(state,seq(count))* max_k(vOld(:)*tr(k,state));
    bestVal = -inf;
    bestPTR = 0;
    if (state == targetState(1)) || (state == targetState(2))
        boost_count(state) = 1;
    else
        boost_count(state) = 0;
    end
    % use a loop to avoid lots of calls to max
    for inner = 1:numStates
      if (count > 200) && (dCount > 20) %&& boost_count(state)
        % add the mean change in probability over last 3 time steps
        vBoost = vStore(inner,count-3:count);
        vBoost = vBoost(~isnan(vBoost));
        vOld_boost = min(0,vOld(inner) + mean(abs(vBoost))); % encourage a positive transition
      else
        vOld_boost = vOld(inner);
      end

      %val = vOld(inner) + logTR(inner,state) ;
      val = vOld_boost + logTR(inner,state);
      if val > bestVal
        bestVal = val;
        bestPTR = inner;
      end
    end
    % save the best transition information for later backtracking
    pTR(state,count) = bestPTR;
    %%% update v - changed from default
    vCalc = 0;
    % calculate for each measure - {distance, orient, vel, accel}
    for meas = 1:size(logE,3) 
      % for every measure using independence assumption
      vCalc = vCalc + logE(state,seq(state,count,meas),meas);      
    end
    v(state) = vCalc + bestVal;
  end
  % store the change in probability
  vStore(:,count) = v-vOld; 
  vOld = v;  
end

% decide which of the final states is post probable
[logP, finalState] = max(v);

% Now back trace through the model
currentState(L) = finalState;
for count = L-1:-1:1
  currentState(count) = pTR(currentState(count+1),count+1);
  if currentState(count) == 0
    error(message('stats:hmmviterbi:ZeroTransitionProbability', currentState( count + 1 )));
  end
end
if customStatenames
  currentState = reshape(statenames(currentState),1,L);
end


