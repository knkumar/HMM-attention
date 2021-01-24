function getOverallStats(model, foldername, fname) 
    % FR - first response
    % NR - number of responses
   

    sNames = dir([foldername,'/',fname]);
    
    result = struct();
    overallStats = struct();
    % 4 time groups - {all, <600, 600-800, >800}
    for grp = 1:4
        %dataMat(grp) = struct();
        dataMat(grp).distMat = [];
        dataMat(grp).l1FRMat = []; dataMat(grp).l1FCMat = []; dataMat(grp).l1NRMat = [];
        dataMat(grp).mdistNRMat = []; dataMat(grp).mdistFRMat = []; dataMat(grp).mdistFCMat = [];
        dataMat(grp).vdistNRMat = []; dataMat(grp).vdistFRMat = []; dataMat(grp).propMat = []; 
        dataMat(grp).propMatc= []; dataMat(grp).propMatCR= []; dataMat(grp).propMatCH= [];
        dataMat(grp).errFRMat = []; dataMat(grp).errFCMat = [];
    end

    targetTypeData = struct('Red',dataMat,'Size',dataMat,'Square',dataMat);
    condData = struct('FF',dataMat,'TF',dataMat,'FT',dataMat,'TT',dataMat);

    fvalsMat = []; saveNames = {};
    transTimes = struct();
    target_map = struct('Red','R','Size','S','Square','Q');

    % calculate overall stats for each subject
    for i = 1:length(sNames)
        sName = sNames(i);
        
        
        sub_name = strrep(regexprep(sName.name,'-center.*',''),'fixed_','');
        sub_name = strrep(sub_name, '-', '_');
        blkTimes = struct('R',[],'S',[],'Q',[]);
        transTimes.(sub_name) = struct('Red',blkTimes,'Square',blkTimes,'Size',blkTimes);

        saveName = sub_name;
        load([sName.folder '/' sName.name]);
        
        if exist('xfin')
            fvalsMat = cat(1,fvalsMat, xfin);
        end
        
        blkNames = fieldnames(summary.respondTargets);
        
        % first correct time matrices for target and foil
        FCmatt = [];
        FCmatf = [];
        
        % get distributon of times from each block
        for blkId = 1:length(blkNames)
            blkName = blkNames{blkId};
            if strcmp(blkName,'agg') | strcmp(blkName, 'Red') | strcmp(blkName, 'Size') | strcmp(blkName, 'Square') | strcmp(blkName, 'aggC')
                continue
            end
            
            targetBlk = summary.respondTargets(1).(blkName);
            foilBlk = summary.respondFoils(1).(blkName);
            
            FCmatt = [FCmatt; targetBlk.firstResC];
            FCmatf = [FCmatf; foilBlk.firstResC];
            
            stim_name = blkName(4:end-1);
            stim = target_map.(stim_name);
            transTimes.(sub_name).(stim_name).(stim) = [transTimes.(sub_name).(stim_name).(stim), targetBlk.firstResC];
            all_stims = {'R','S','Q'};

            for letter_id = 1:length(all_stims)
                letter = all_stims{letter_id};
                if sum(letter == foilBlk.stim)
                    stim = letter;
                    transTimes.(sub_name).(stim_name).(stim) = [transTimes.(sub_name).(stim_name).(stim), foilBlk.firstResC(foilBlk.stim==stim)];
                end
            end

        end
        for timeGrp = 1:4
            % disp(['Time: ',num2str(timeGrp)])
            targetData = summary.respondTargets(timeGrp).agg;
            foilData = summary.respondFoils(timeGrp).agg;
            
            
            [targetData, targetNan] = makeChecks(targetData);
            [foilData, foilNan] = makeChecks(foilData);
            datanan = targetNan | foilNan;
            
            blocks = targetData.blocks;
            % get only the first block
            oneBlk = cellfun(@(x) ~isempty(x), strfind(blocks,'one'));
            saveNames{i} = saveName;
            % disp(['Block:',blkName,'-fillData'])
            [result, dataMat(timeGrp)] = fillData(targetData, foilData, saveName, result, dataMat(timeGrp), FCmatt, FCmatf, oneBlk);
            
            % NK : looking at specific target types
            targetTypes = {'Red','Size','Square'};
            for ttype = 1:length(targetTypes)
                % disp(targetTypes{ttype})
                targetName = targetTypes{ttype};
                TargetData = summary.respondTargets(timeGrp).(targetName);
                FoilData = summary.respondFoils(timeGrp).(targetName);
                [ignore, targetTypeData.(targetName)(timeGrp)] = fillData(TargetData, FoilData, saveName, [], targetTypeData.(targetName)(timeGrp), FCmatt, FCmatf, oneBlk);
            end
            
            conditional = {'FF','FT','TF','TT'};
            for ttype = 1:length(conditional)
                % disp(conditional{ttype})
                cond_name = conditional{ttype};
                target_conditional = summary.respondTargets(timeGrp).aggC.(cond_name);
                foil_conditional = summary.respondFoils(timeGrp).aggC.(cond_name);
                [ignore, condData.(cond_name)(timeGrp)] = fillData(target_conditional, foil_conditional, saveName,[], condData.(cond_name)(timeGrp), FCmatt, FCmatf, oneBlk);
            end
        end
    end
    save([foldername,'/Overall',regexprep(model,'*',''),'.mat'],'saveNames', 'result', 'fvalsMat', 'dataMat', 'targetTypeData', 'condData', 'transTimes');
end

function [result,dataMat] = fillData(targetData, foilData, saveName, result, dataMat, FCmatt, FCmatf, oneBlk)
    % oneBlk - get only the first block
    % get means from target and foil data
    
    distFRt = [targetData(:).meanFR];
    distFRf = [foilData(:).meanFR];

    distFCt = [targetData(:).meanFC];
    distFCf = [foilData(:).meanFC];

    distNRt = [targetData(:).meanNR];
    distNRf = [foilData(:).meanNR];

    N = length(distFRt);

    tFRFit = glmfit(1:N,distFRt);
    fFRFit = glmfit(1:N,distFRf);
    % fit a line through target number of transitions
    tNRFit = glmfit(1:N,distNRt);
    fNRFit = glmfit(1:N,distNRf);

    if ~isempty(result)
        result.(saveName).dist1 = distFRt-distFRf;
        result.(saveName).dist2 = distNRt-distNRf;
        result.(saveName).mdist = mean(sqrt((distFRt-distFRf).^2+(distNRt-distNRf).^2));
        result.(saveName).tFRFit = tFRFit;
        result.(saveName).fFRFit = fFRFit;
        result.(saveName).tNRFit = tNRFit;
        result.(saveName).fNRFit = fNRFit;
        result.(saveName).l1FR = (tFRFit(2) - fFRFit(2));
        result.(saveName).l1NR = (tNRFit(2) - fNRFit(2));
        result.(saveName).FCmatt = FCmatt;
        result.(saveName).FCmatf = FCmatf;
        % disp(saveName{1})
        % size(mean(targetData(:).meanpropC))
    else
        oneBlk = ones(size(distFRt));
    end


    dataMat.distMat = cat(2, dataMat.distMat, mean(sqrt((distFRt-distFRf).^2+(distNRt-distNRf).^2)) );
    % dprime for any transition
    dataMat.l1FRMat = cat(2, dataMat.l1FRMat, (tFRFit(2) - fFRFit(2)) ) ;

    % NK: 1/15/2018 - added a dprime in time for correct transitions
    dataMat.l1FCMat = cat(2, dataMat.l1FCMat, mean(distFCt-distFCf));

    dataMat.l1NRMat = cat(2, dataMat.l1NRMat, (tNRFit(2) - fNRFit(2)));
    % disp('Size of num FC')
    % disp([targetData(:).num_FC])
    dataMat.mdistFCMat = cat(2, dataMat.mdistFCMat, [mean(distFCt); mean(distFCf); ...
                            mean(distFCt(oneBlk)); mean(distFCf(oneBlk)); mean(distFCt(~oneBlk)); mean(distFCf(~oneBlk)); ...
                            sum([targetData(:).num_FC]); sum([foilData(:).num_FC]) ]);
    dataMat.mdistFRMat = cat(2, dataMat.mdistFRMat, [mean(distFRt); mean(distFRf); ...
                            mean(distFRt(oneBlk)); mean(distFRf(oneBlk)); mean(distFRt(~oneBlk)); mean(distFRf(~oneBlk)) ]);
    dataMat.mdistNRMat = cat(2, dataMat.mdistNRMat, [mean(distNRt); mean(distNRf); ...
                            mean(distNRt(oneBlk)); mean(distNRf(oneBlk)); mean(distNRt(~oneBlk)); mean(distNRf(~oneBlk)) ]);
    dataMat.vdistFRMat = cat(2, dataMat.vdistFRMat, [var(distFRt); var(distFRf)]);
    dataMat.vdistNRMat = cat(2, dataMat.vdistNRMat, [var(distNRt); var(distNRf)]);

    errFRt = [targetData(:).errFR];
    errFRf = [foilData(:).errFR];
    dataMat.errFRMat = cat(2, dataMat.errFRMat, [ mean(errFRt); mean(errFRf) ]);

    errFCt = [targetData(:).errFC];
    errFCf = [foilData(:).errFC];
    dataMat.errFCMat = cat(2, dataMat.errFCMat, [ mean(errFCt); mean(errFCf) ]);
    % disp('Size of num prop')
    % disp(targetData(:).num_prop)
    dataMat.propMat = cat(2, dataMat.propMat, [mean([targetData(:).meanprop]); mean([foilData(:).meanprop]); mean([targetData(:).errNR]);...
                             mean([foilData(:).errNR]); sum([targetData(:).num_prop]); sum([foilData(:).num_prop]) ]);
    dataMat.propMatc = cat(2, dataMat.propMatc, [mean([targetData(:).meanpropC]); mean([foilData(:).meanpropC]); mean([targetData(:).errpropC]);...
                             mean([foilData(:).errpropC]); sum([targetData(:).num_propC]); sum([foilData(:).num_propC]) ]);
    dataMat.propMatCR = cat(2, dataMat.propMatCR, [mean([targetData(:).meanpropCR]); mean([foilData(:).meanpropCR]); mean([targetData(:).errpropCR]); mean([foilData(:).errpropCR]) ]);
    dataMat.propMatCH = cat(2, dataMat.propMatCH, [mean([targetData(:).meanpropCH]); mean([foilData(:).meanpropCH]); mean([targetData(:).errpropCH]); mean([foilData(:).errpropCH]) ]);

end


function [data,dataNan] = makeChecks(data)
    % check data for nans and modify to zero
    % input  : data with nans
    % output : data without nans
    %          location of nans
    data.meanFR(isnan(data.meanFR)) = 0;
    data.meanFC(isnan(data.meanFC)) = 0;
    data.meanpropC(isnan(data.meanpropC)) = 0;
    data.meanprop(isnan(data.meanprop)) = 0;
    data.meanRes(isnan(data.meanRes)) = 0;
    dataNan = isnan(data.meanFR);
end
