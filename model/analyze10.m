function [targets,hand,sim,hmm] = analyze10(targets,hand,sim,hmm)

    %% HMM properties
    if ~exist('hmm','var')
        x = [.7, .1, ...
            0, 1, 1, 1, ... % dist
            1, 1, ... % orient (1, 1, 0, 1,)
            1, 1, 0.2, ... % vel
            1, 1]; % acc
        hmm = getHMM_null(x);
    end
    % add center state
    targets.tLocs(4,:) = fliplr(targets.Resolution/2);

    generate_emiss

        function generate_emiss
            names = hmm.emit.names;
            % compute emission prob for each variable
            for n = 1:length(names) % {distance, orient, vel, accel}
                hmm.emit.(names{n}).prob = [];
                % run for each state in model
                for st = 1:length(hmm.trans)
                    xs = linspace(hmm.emit.(names{n}).range(1),hmm.emit.(names{n}).range(2),hmm.n.bins);
                    dist = hmm.emit.(names{n}).dist(st,:);
                    % generate 2-parameter pdf for xs values
                    ps = pdf(dist{1},xs,dist{2},dist{3});
                    ps(isinf(ps)) = 1;
                    % normalize the density to sum to 1
                    ps = ps ./ sum(ps);
                    ps(isnan(ps)) = 0;
                    psArrange = ps([hmm.emit.(names{n}).offset:length(ps),1:hmm.emit.(names{n}).offset-1]);
                    if (strcmp(names{n},'distance'))
                        ps(1) = ps(2);
                    end
                    hmm.emit.(names{n}).xs = xs;
                    hmm.emit.(names{n}).prob = [hmm.emit.(names{n}).prob;psArrange];
                end
                % store all emission probabilities
                hmm.emit.all.prob = cat(3,hmm.emit.all.prob,hmm.emit.(names{n}).prob);
                hmm.emit.all.prob(hmm.emit.all.prob==0) = 1e-20; %add infinitesmal probability
            end
        end



    %% get probabilities from data
    calc_hand
    calc_obj
    calc_emiss

    % Calculate velocity and acceleration
        function calc_hand
            hand.calc.x = hand.data.x';
            hand.calc.y = hand.data.y';
            hand.calc.dx = diff([0;hand.calc.x]);
            hand.calc.dy = diff([0;hand.calc.y]);
            hand.calc.velocity = sqrt(hand.calc.dx.^2 + hand.calc.dy.^2);
            hand.calc.velocity_l1 = abs(hand.calc.dx) + abs(hand.calc.dy);
            hand.calc.velocity_lmax = max(abs(hand.calc.dx), abs(hand.calc.dy));

            hand.calc.acceleration = diff([0;hand.calc.velocity]);
            hand.calc.angle = atan2d(hand.calc.dy,hand.calc.dx);

            hand.calc.dxNoise = hand.calc.dx +.001*randn(size(hand.calc.dx));
            hand.calc.dyNoise = hand.calc.dy +.001*randn(size(hand.calc.dy));
        end

    % calculate target distance and orientation
        function calc_obj
            targets.calc.position.x = ones(length(hand.calc.dxNoise),1) * targets.tLocs(:,1)'; % target location x vector
            targets.calc.position.y = ones(length(hand.calc.dyNoise),1) * targets.tLocs(:,2)'; % target location y vector

            targets.calc.dx = targets.calc.position.x - (hand.calc.x * ones(1,floor(hmm.n.states/2)));
            targets.calc.dy = targets.calc.position.y - (hand.calc.y * ones(1,floor(hmm.n.states/2)));

            targets.calc.distance = sqrt(targets.calc.dx.^2 + targets.calc.dy.^2);

            targets.calc.distance(:,4) = targets.calc.distance(:,4)*2;

            targets.calc.angles = nan(length(hand.calc.x),ceil(hmm.n.states/2));
            hVec = [hand.calc.dx,hand.calc.dy]';
            for T = 1:floor(hmm.n.states/2)
                oVec = [targets.calc.dx(:,T),targets.calc.dy(:,T)]';
                targets.calc.angles(:,T) = acosd(dot(hVec,oVec)./(sqrt(dot(hVec,hVec)).*sqrt(dot(oVec,oVec))));
            end
            targets.calc.angles(isnan(targets.calc.angles)) = 0;
            targets.calc.position.current = [targets.tLocs(targets.data.state.current,1),...
                targets.tLocs(targets.data.state.current,2)]';
            tNan = find(isnan(targets.data.state.last));
            targets.data.state.last(tNan) = 1;
            targets.calc.position.last= [targets.tLocs(targets.data.state.last,1),...
                targets.tLocs(targets.data.state.last,2)]';
            targets.calc.position.last(:,tNan) = nan;
            targets.data.state.last(tNan) = nan;
        end

    % Calculate emission probablity
        function calc_emiss
            hmm.calc.velocity = ...
                max(1,min(hmm.n.bins,round(hand.calc.velocity'/(diff(hmm.emit.velocity.range)/hmm.n.bins))));
            hmm.calc.acceleration = ...
                max(1,min(hmm.n.bins,round([hand.calc.acceleration]'/(diff(hmm.emit.acceleration.range)/hmm.n.bins))));
            nT = size(targets.tLocs,1);
            for T = 1:nT
                hmm.calc.distance(T,:) = ...
                    max(1,min(hmm.n.bins,round(targets.calc.distance(:,T)'/(diff(hmm.emit.distance.range)/hmm.n.bins))));
                hmm.calc.orientation(T,:) = ...
                    max(1,min(hmm.n.bins,round(targets.calc.angles(:,T)'/(diff(hmm.emit.orientation.range)/hmm.n.bins))));
            end
            hmm.calc.orientation(:,hmm.calc.velocity==1) = 1;

            names = {'distance','orientation','velocity','acceleration'};
            hmm.calc.all = [];
            for n = 1:length(names)
                hmm.calc.all = cat(3,hmm.calc.all,hmm.calc.(names{n})([hmm.emit.(names{n}).order],:));
            end
        end

    [hmm.est.state,hmm.est.vStore] = hmmviterbi_multi2(hmm.calc.all,hmm.trans,hmm.emit.all.prob, hand.data.trialID, hand.data.trial);

    hmm.est.objState = nan(size(hmm.est.state));
    hmm.est.objState(hmm.est.state == 1 | hmm.est.state == 2) = 1;
    hmm.est.objState(hmm.est.state == 3 | hmm.est.state == 4) = 2;
    hmm.est.objState(hmm.est.state == 5 | hmm.est.state == 6) = 3;
    hmm.est.objState(hmm.est.state == 7) = 0;
    hmm.est.objState(hmm.est.state == 8) = 4;

end
