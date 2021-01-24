% function to plot hmm model predictions
function plotModel(targets,hmm, save_name, plot_resolution, plot_title, save_dir)
    periods = targets.data.period;
    block = targets.data.block;
    
    for i = 1:length(block)
        blocks{i} = [periods{i} '_' block{i}];
    end
    
    xres = plot_resolution(1);
    yres = plot_resolution(2);
    
    [~,idx] = unique(blocks);
    uniqueBlocks = blocks(sort(idx));
    % plot each block
    for i = 1:length(uniqueBlocks)
        
        fig     = figure('Visible','Off');
        hold on;
        fig.PaperPositionMode = 'manual';
        fig.PaperType = '<custom>';
        fig.PaperUnits  = 'points';
        fig.PaperSize   = [xres yres];
        fig.Position    = [0 0 xres yres];
        fig.OuterPosition = [0 0 xres yres];
        currentBlk      = uniqueBlocks{i};
        % Identify the indexes in data for current Block
        blkIds          = strcmp(blocks,currentBlk); 
        actualStates    = targets.data.state.current(blkIds);
        % Generate unique trials
        trials          = targets.data.trial(blkIds);
        trialNum        = 1:length(trials);
        [~,idx]         = unique(trials);
        uniqueTrials    = trials(sort(idx));
        estStates = hmm.est.objState(blkIds);
        
        n = ceil(length(trials)/2000);
        oldk=0;
        % plot trials
        for j = 1:length(uniqueTrials)
            currTrial   = uniqueTrials{j};
            trialIds    = strcmp(trials,currTrial);
            currTrialNum = trialNum(trialIds);
            newk = ceil(max(currTrialNum)/2000);
            % create a new subplot every 1000 time steps
            if oldk~=newk
                oldk=newk;
                r = ceil(n/2);
                ax(oldk) = subplot(n,1,oldk);
                xbar = min(currTrialNum):oldk*2000;
                patch_plot(xbar,ax(oldk));
            end
            x=[currTrialNum fliplr(currTrialNum)];
            y1=actualStates(currTrialNum);
            y2=estStates(currTrialNum);
            if isletter(currTrial(2))
                p1 = patch('Faces',1:length(x), 'Vertices',[x' ,[y1-0.3, y1+0.3]'],'FaceColor','red','EdgeColor','black','LineWidth',0.5);
            else
                p2 = patch('Faces',1:length(x), 'Vertices',[x' ,[y1-0.3, y1+0.3]'],'FaceColor','blue','EdgeColor','black','LineWidth',0.5);
            end
            y2(end)=nan;
            p3 = patch(currTrialNum,y2,'black','EdgeColor','black','FaceColor','None','Marker','None','LineWidth',0.5);
            ylim([-1 5]);
        end
        
        %% Plot Legend
        h=legend([p1,p2,p3],'Foils','Targets','Estimated States');
        set(gcf,'Units','normalized');
        rect=[0.75, 0.94, .1, .025]; set(h,'Position',rect,'FontSize',10,'LineWidth',0.1);
        
        ha = axes('Position',[0 0 1 0.98],'Xlim',[0 1],'Ylim',[0 
                 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.65, 1,strcat('\bf',plot_title),'HorizontalAlignment', ...
                'right','VerticalAlignment', 'top','Color','black','FontSize',12); 
        text(0.5, 0.05 ,'\bf Time in 10 ms intervals','HorizontalAlignment', ...
                'center','VerticalAlignment', 'top','Color','black','FontSize',10);   
            
        [im_hatch,colorlist] = applyhatch_pluscolor(fig,'x',1,[],[255 0 0],'300',1,1); % apply patch for red color
        [im_hatch,colorlist] = applyhatch_pluscolor(im_hatch,'w',1,[],[0 0 255],'300',3,1); % apply patch for blue color
        iptsetpref('ImshowAxesVisible','off');

        if ~exist(save_dir, 'dir')
           mkdir(save_dir) 
        end
        
        fName = fullfile(save_dir, strcat(save_name,currentBlk));
        imwrite(rgb2gray(im_hatch),strcat(fName,'.png'),'png','XResolution',xres,'YResolution',yres);
        
        close(fig); close all;
    end
    
end

function patch_plot(xbar,ax)
    xbars = [xbar fliplr(xbar)];
    y=zeros(1,length(xbar));
    
    ax.TickLength = [0.005 0.005];
    ax.XTickMode='manual'; ax.YTickMode;'manual';
    ax.YTick=[0 1 2 3 4]; 
    ax.YTickLabel={'','State 1','State 2','State 3',''};
    ax.XLim=[min(xbar) ceil(max(xbar)/1000)*1000];
    ax.XMinorTick='on'; ax.XTickLabelRotation=30; ax.TickDir='in';
    ax.XTick=min(xbar):50:ceil(max(xbar)/1000)*1000; ax.Box='off';
    ax.FontSize=6;
    ax.Position(1) = ax.Position(1) * 0.75;
    ax.Position(3) = ax.Position(3) * 1.1;
    ax.Position(4) = ax.Position(4) * 0.9;
    ax.LineWidth = 0.1;
end