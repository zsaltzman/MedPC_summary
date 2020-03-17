%% rick2019_12_iterative_v2
%built from rick2019_11_iterative


%v2
    %1 sec bins for latencies and changed graphing to autoupdate to bins
    %latencybin(2:end)

%Updating notes:
%changed the importing to have the folder and file be changed in
%the beginning

%changed the way the columns are indexed from raw by using string
%compare with the header, to avoid possible slipups

%%

%To Do
%update input file
%update file output names
%update row/column names
%eventually update the way I read in timestamps to not be dependent on
%changing the number of columns to pull and instead read to column
%heading to determine that (using strcmp and matching the first two
%characeters)

clear; clc

%mice in each group
%more mice are in groups but left out shitters
saline_mice = [5 12	16 19 26 37];
mec_mice = [6 9 13 20 24 27 38];
scop_mice = [3	7 10 14	17 21 28 32	35 39];
mec_scop_mice = [4 8 11 15 18 22 25 29 36 40];

groupkey = {'Saline', 'Mec', 'Scop', 'Mec+Scop'};

salinedataindex = 0;
mecdataindex = 0;
scopdataindex = 0;
mec_scopdataindex = 0;

groupdata = {'salinedata'; 'mecdata'; 'scopdata'; 'mec_scopdata'};
groupdataindex = {'salinedataindex', 'mecdataindex', 'scopdataindex', 'mec_scopdataindex'};


subsetday = 20; %last day to run for analyzing subset of days


salinemouseorder = [];
mecmouseorder = [];
scopmouseorder = [];
mec_scopmouseorder = [];

folder = 'C:\Users\User\Google Drive\2019-12 Antagonist Pilot\' ;
inputfile = '2019-12 MATLAB Full';

codename = 'rick2019_12_iterative';
outputfolder = 'C:\Users\User\Google Drive\2019-12 Antagonist Pilot\2019-12 MATLAB\';
graphfolder = 'C:\Users\User\Google Drive\2019-12 Antagonist Pilot\2019-12 MATLAB\Graphs\';

outputlatfile = '2019-12 MATLAB Mean Lat, Rew Lat, Acq Time.xlsx';
outputmatfile = '2019-12 MATLAB Variables.mat';


%set the variable letters that you're pulling
%Correct = B, Inactive = D, Receptacle = G, Reward = H, Tone on = K
%Tone off = L, Incorrect = R, Intervals used = S
variable_letters = ["B(" , "D(" , "G(" , "H(" , "K(" , "L(" , "R(" , "S("];

%graphing/binning parameters to change
latencyedges = 0:1:10;

%% Import the data
[~, ~, raw] = xlsread([folder inputfile]);
%only imported Cued and CuedTO

%cut off column headings
%but keep the summary column headings separate
medheadsum = raw(1,1:17);
medheader = raw(1,:);
raw = raw(2:end,:);




%Identify all unique mouse ID numbers and save in variable
raw_mouse_ID = unique(cell2mat(raw([1:end],1)));


%trim to just eyfp (that ran against chr2, not the arch eyfp control) and chr2 mice
all_mouse_idx = 0;

for mouse_ID_idx = 1:size(raw_mouse_ID,1)
    %remove this if statement (but leave the two lines in between), if I want to include the excluded mice) 
    if any(raw_mouse_ID(mouse_ID_idx) == saline_mice) || any(raw_mouse_ID(mouse_ID_idx) == mec_mice) || any(raw_mouse_ID(mouse_ID_idx) == scop_mice)|| any(raw_mouse_ID(mouse_ID_idx) == mec_scop_mice)
        all_mouse_idx = all_mouse_idx +1;
        all_mouse_ID(all_mouse_idx,1) = raw_mouse_ID(mouse_ID_idx);
    end
end

%Create output variables for aggregate mean latency and acquisition time
%Number of days + 1 (for mouse # headers) , number of mice
%May cause errors if all mice don't run the same number of days
AllMeanLat = cell(round((size(raw,1)/size(all_mouse_ID,1)+1)),size(all_mouse_ID,1));
AllRewMeanLat = cell(round((size(raw,1)/size(all_mouse_ID,1)+1)),size(all_mouse_ID,1));
AllAcqTime = cell(round((size(raw,1)/size(all_mouse_ID,1)+1)),size(all_mouse_ID,1));

%cycle through each mouse
for num = 1:length(all_mouse_ID)
    
    %Define mouse_ID number for the run of the for loop
    mouse_ID = all_mouse_ID(num);
    
    %Cut down raw to just current mouse
    raw_mouse = raw(cell2mat(raw(:,1)) == mouse_ID,:);
    
    
       %select group mouse belongs to
        if strcmp(raw_mouse{1,6},'Saline')
            group = 1;
            salinedataindex = salinedataindex + 1;
            salinemouseorder(salinedataindex) = mouse_ID;
                         
        elseif strcmp(raw_mouse{1,6},'Mecamylamine')
            group = 2; 
            mecdataindex = mecdataindex + 1;
            mecmouseorder(mecdataindex) = mouse_ID;
            
        elseif strcmp(raw_mouse{1,6},'Scopolamine')
            group = 3; 
            scopdataindex = scopdataindex + 1;
            scopmouseorder(scopdataindex) = mouse_ID;
            
        elseif strcmp(raw_mouse{1,6},'Mec+Scop')
            group = 4; 
            mec_scopdataindex = mec_scopdataindex + 1;
            mec_scopmouseorder(mec_scopdataindex) = mouse_ID;
            
            
        end
    
    %% cycle through each day
    
    for row = 1:size(raw_mouse,1)
        % if you want to do a subset of days
        %         for row = 1:subsetday
        
        
        for variable=1:8
            %% Pull the variables from raw into their own rows of cell raw1
            
            %Find all the days that are the given variable
            array_header = strfind(medheader,variable_letters{variable});
            
            %find the first non-empty cell and set that index to TO
            first = find(~cellfun(@isempty,array_header),1);
            last = find(~cellfun(@isempty,array_header),1, 'last');
            
            raw1 = raw_mouse(row,first:last);
            
            
            raw1(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw1)) = {''};
            
            %% Replace non-numeric cells with NaN
            R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw1); % Find non-numeric cells
            raw1(R) = {NaN}; % Replace non-numeric cells
            
            %% Create output variable
            data{row,variable}= reshape([raw1{:}],size(raw1));
            
            
            %% Clear temporary variables
            clearvars raw1 R ;
            
        end
        
        %% Summary values
        %Note: will be adding to data cell after timeouts to keep
        %consistent with future codes (though, this is in a different place
        %than the FP stuff, which has medsum in the beginning)
        medsumdata=raw_mouse(row,1:17);
        medsum{row,1} = [medheadsum;medsumdata];
        
        
    end
    clear raw_mouse
    
    
    %% Calculate Time of Timeouts (search through improper np's that triggered a timeout i.e. happened 5 sec after last one)
    % ! won't be needed for future versions that have this built in as a
    %   variable
    
    timeout = cell(size(data,1),1);
    
    for datarow = 1:size(data,1)
        timeout{datarow,1}(1) = data{datarow,7}(1);
        timeoutidx = 1;
        
        for incorrect = 2:size(data{datarow,7},2)
            if data{datarow,7}(incorrect) - timeout{datarow,1}(timeoutidx) >= 5
                timeoutidx = timeoutidx + 1;
                timeout{datarow,1}(timeoutidx) = data{datarow,7}(incorrect);
                
            end
            
        end
        
    end
    
    
    %append timeout to data
    data = [data timeout];
    
    
    %% Omissions / tone miss
    omission = cell(size(data,1),1);
    
    for datarow = 1:size(data,1)
        omissionidx = 0;
        
        for tone = 1:size(data{datarow,5},2)
            if data{datarow,6}(tone) - data{datarow,5}(tone) == 10
                omissionidx = omissionidx + 1;
                omission{datarow,1}(omissionidx) = data{datarow,5}(tone);
                
            end
            
        end
        
    end
    
    
    %append timeout to data
    data = [data omission];
    
    
    
    %% Add medsum to data cell
    data = [data medsum];
    
    %% Calculate latency (cue on to Reward(First Proper NP of trial), col 9) and training day mean latency (col 10), and concat to data. Row = training day
    
    %Changed ProperANP to rewards from previous version because ProperANPs now
    %colllect any extra Proper ANPs made (within the 2 sec the tone can still play after the first), not just first one
    latency = cell(size(data,1),4);
    for datarow=1:size(data,1)
        
        for Reward=1:size(data{datarow,4},2)
            if data{datarow,4}(Reward)>0
                toneidx = find(data{datarow,4}(Reward)>data{datarow,5},1, 'last');
                if data{datarow,4}(Reward) - data{datarow,5}(toneidx) <= 10
                    latency{datarow,1}(Reward) = data{datarow,4}(Reward)- data{datarow,5}(toneidx);
                end
            end
        end
        
        
        
        %calc average latency
        latency{datarow,2} = cellfun(@mean,latency(datarow));
        
        %break into bins
        
        latency{datarow,3} = histcounts(latency{datarow,1},latencyedges);
        
        %normalize counts
        latency{datarow,4} = latency{datarow,3}(1,:)/sum(latency{datarow,3}(1,:));
        
    end
    
    %add latency to data cell
    data = [data latency];
    
    %add mean latency to aggregate
    AllMeanLat{1,num} = mouse_ID;
    if size(AllMeanLat,1) == size(latency,1)+1
        AllMeanLat(2:end,num) = latency(:,2);
        
    else
        latency{size(AllMeanLat,1)-1,2} = NaN;
        AllMeanLat(2:end,num) = latency(:,2);
    end
    
    
    
    clear datarow ProperANP toneid Reward
    
    %% Calculate rewlatency (Reward to Receptacle latency col 12) and training day mean latency (col 13), and concat to data. Row = training day
    %not finished writing this part
    
    rewlatency = cell(size(data,1),4);
    for datarow=1:size(data,1)
        
        
        for Reward=1:size(data{datarow,4},2)
            if data{datarow,4}(Reward)>0
                recepidx = find(data{datarow,4}(Reward)<data{datarow,3},1);
                if data{datarow,3}(recepidx)-data{datarow,4}(Reward)<10
                    rewlatency{datarow,1}(Reward) = data{datarow,3}(recepidx) - data{datarow,4}(Reward);
                end
            end
        end
        
        %calculate mean of all rewlatencies
        rewlatency{datarow,2} = cellfun(@mean,rewlatency(datarow));
        
        %break into bins
        
        rewlatency{datarow,3} = histcounts(rewlatency{datarow,1},latencyedges);
        
        %normalize counts
        rewlatency{datarow,4} = rewlatency{datarow,3}(1,:)/sum(rewlatency{datarow,3}(1,:));
    end
    
    %add latency to data cell
    data = [data rewlatency];
    
    %add mean latency to aggregate
    AllRewMeanLat{1,num} = mouse_ID;
    if size(AllRewMeanLat,1) == size(rewlatency,1)+1
        AllRewMeanLat(2:end,num) = rewlatency(:,2);
        
    else
        rewlatency{size(AllRewMeanLat,1)-1,2} = NaN;
        AllRewMeanLat(2:end,num) = rewlatency(:,2);
    end
    
    clear datarow ProperANP toneid  Reward
    %% Time to reach 30 rewards, acqtime (col 11)
    Acqtime = cell(size(data,1),1);
    
    for datarow = 1:size(data,1)
        if isnan(data{datarow,4}(30)) || data{datarow,4}(30) == 0
            Acqtime{datarow,1} = 1800;
            
        else
            Acqtime{datarow,1} = data{datarow,4}(30);
            
        end
    end
    
    %add Acqtime to data cell
    data = [data Acqtime];
    
    %add Acqtime to aggregate
    AllAcqTime{1,num} = mouse_ID;
    
    if size(AllAcqTime,1) == size(Acqtime,1)+1
        AllAcqTime(2:end,num) = Acqtime;
        
    else
        Acqtime{size(AllAcqTime,1)-1,1} = NaN;
        AllAcqTime(2:end,num) = Acqtime;
    end
    
    clear datarow
    
    %% Replace cells that are all 0's with same size cells of NaNs
    %this is to not get wonky stuff happening when I graph cmd
    %plots
    
    for datarow = 1:size(data,1)
        for datacol = 1:size(data,2)
            %skip datacol 11 becuase it's medsum
            if datacol == 11
                continue
                
            elseif sum(data{datarow,datacol}) == 0
                data{datarow,datacol} = NaN(size(data{datarow,datacol},1), size(data{datarow,datacol},2));
            end
        end
    end
    
    clear datarow datacol
    
    
    %% Action Hist count
    %do a histcount for variables except tone offa intervals used
    %these will have actionedges = 0:300:1800;
    
    %     %previously did days as individual cell rows.
    %     actioncounts = cell(size(data,1),10);
    %     %normalized, i.e. percent of total counts for that session
    %     normactioncounts = cell(size(data,1),10);
    %     actionedges = 0:300:1800;
    %
    %     for datarow = 1:size(data,1)
    %         %skip datacol 6 (tone off) and 8 (intervals used)
    %
    %         for datacol = 1:10
    %             if  datacol ==6 || datacol==8
    %                 continue
    %
    %             else
    %                 actioncounts{datarow,datacol} = histcounts(data{datarow,datacol},actionedges);
    %                 normactioncounts{datarow,datacol} = actioncounts{datarow,datacol}(1,:)/sum(actioncounts{datarow,datacol}(1,:));
    %             end
    %         end
    %     end
    
    
    actioncounts = cell(1,10);
    %normalized, i.e. percent of total counts for that session
    normactioncounts = cell(1,10);
    actionedges = 0:300:1800;
    
    for datarow = 1:size(data,1)
        %skip datacol 6 (tone off) and 8 (intervals used)
        
        for datacol = 1:10
            if  datacol ==6 || datacol==8
                actioncounts{1,datacol} = NaN(size(data,1),6);
                normactioncounts{1,datacol} = NaN(size(data,1),6);
                
                continue
                
            else
                actioncounts{1,datacol}(datarow,:) = histcounts(data{datarow,datacol},actionedges);
                normactioncounts{1,datacol}(datarow,:) = actioncounts{1,datacol}(datarow,:)/sum(actioncounts{1,datacol}(datarow,:));
            end
        end
    end
    
    
    %% Add data to group cells
    %Note: before data row/column labels
    %standard adding of mouse_ID and data
    eval([groupdata{group},'{1,' [groupdata{group} 'index'] '} = mouse_ID;' ]);
    eval([groupdata{group},'{2,' [groupdata{group} 'index'] '} = data;']);
    
    %add actioncount (the histcounts version of each action),
    %normactioncount, and histcount(latencies) (both count and norm)
    eval([groupdata{group},'{3,' [groupdata{group} 'index'] '} = actioncounts;']);
    eval([groupdata{group},'{4,' [groupdata{group} 'index'] '} = normactioncounts;']);
    eval([groupdata{group},'{5,' [groupdata{group} 'index'] '} = cell2mat(data(:,14));']);
    eval([groupdata{group},'{6,' [groupdata{group} 'index'] '} = cell2mat(data(:,15));']);
    eval([groupdata{group},'{7,' [groupdata{group} 'index'] '} = cell2mat(data(:,18));']);
    eval([groupdata{group},'{8,' [groupdata{group} 'index'] '} = cell2mat(data(:,19));']);
    
    %add group row labels before saving
    
    
    %% Add row/column labels
    
    %All days ro labels
%     rowlabel = {'Cued 1'  'Cued 2'  'Cued 3'  'Cued 4'  'Cued TO 1' 'Cued TO 2' 'Cued TO 3'  'Cued TO 4'  'Cued TO 5'  'Cued TO 6' 'Cued TO 7' ...
%         'Cued TO 8' 'Cued TO 9' 'Cued TO 10' 'Cued TO 11' 'Cued TO 12' 'Ext 1' 'Ext 2' 'Ext 3'}';
    
        rowlabel = {'Cued 1' ; 'Cued 2'; 'Cued 3'; 'Cued 4' ...
    ; 'Cued TO 1'; 'Cued TO 2'; 'Cued TO 3'; 'Cued TO 4';'Cued TO 5'; 'Cued TO 6';'Cued TO 7';'Cued TO 8'...
    ;'Cued TO 9';'Cued TO 10'; 'Cued TO 11'; 'Cued TO 12'; 'Cued TO Switch 13'; 'Ext 1'; 'Ext 2'; 'Ext 3';...
    'Switcheroo 1'; 'Switcheroo 2'; 'Switcheroo 3'; 'Switcheroo 4'};
    
    data = [rowlabel(1:size(data,1)) data];
    
    columnlabel = {strcat('mouse',num2str(mouse_ID)) 'Correct' 'Inactive' 'Receptacle' 'Reward' ...
        'Tone On' 'Tone Off' 'Incorrect' 'Intervals used (not a timestamp)' 'Timeout' 'Omission' 'Med Summary' 'Tone-Correct Latency' 'Average Tone-Correct Latency' ...
        'Tone-Correct Hist Count' 'Tone-Correct Norm Hist Count' 'Reward-Rec Latency' 'Avaerage Reward-Rec Latency' 'Rew-Rec Hist Count' 'Rew-Rec Norm Hist Count' 'Acquisition Time (time to 30 rewards)'};
    
    data = [columnlabel ; data];
    
    
    clear rowlabel
    %% Create unique variable ID for data and clear variables
    
    eval(['data_', num2str(mouse_ID), '= data;']);
    
    
    
    
    clear datarow Acqtime data dataname variable row
    
    
    
    
    
end


%% Calculate average for groups

for group = 1:size(groupdata,1)
    
    eval(['tempavgdata = ', groupdata{group}, ';']);
    
    Rick_Cat_Code;
    
    
    eval([ groupdata{group}, ' = [', groupdata{group}, ',avgdata];']);
    
end

clear tempavgdata avgdata
%% Plot

for group = 1:size(groupdata,1)
    
    eval(['tempplotdata = ', groupdata{group}, ';']);
    
    for mouse = 1:size(tempplotdata,2)
        for row = [3 5 7] % do by count/norm pair of rows
            
            if row == 3
                
                for action = 1:size(tempplotdata{row,mouse},2)
                    
                    if action == 6 || action == 8
                        %skip the non-actions
                        continue
                        
                    else
                        
                        actionbin = 5:5:30;
                        counts = cell2mat(tempplotdata{row,mouse}(:,action));
                        normcounts = cell2mat(tempplotdata{row+1,mouse}(:,action));
                        
                        
                        figure
                        % tlt=tiledlayout(1,2,'TileSpacing','compact');
                        tlt=tiledlayout(1,2);
                        nexttile(2)
                        imagesc(actionbin, 1, normcounts, [0,1]); % this is the heatmap
                        set(gca,'YDir','normal') % put the trial numbers in better order on y-axis
                        %fix title
                        %title([variables{action} ' Heat Map'],'fontsize',16)
                        colormap(gray)
                        ylabel('Training Day','fontsize',14)
                        xlabel('Time Bin (min)','fontsize',14)
                        cb = colorbar;
                        % !Doesn't scale in window properly, fix later
                        % ylabel(cb, 'Probability','fontsize',16)
                        xticks([0:5:30]);
                        
                        
                        nexttile(1)
                        bstack = bar(counts, 'stacked','FaceColor', 'flat');
                        % barcolor = colormap(parula(6));
                        bincolor = parula(size(counts,2));
                        for binnum = 1:size(counts,2)
                            bstack(binnum).CData = bincolor(binnum,:);
                        end
                        xlabel('Training Day','fontsize',14);
                        ylabel('Raw Count','fontsize',14);
                        
                        
                        if ~ischar(tempplotdata{1,mouse})
                            temptitle = [num2str(tempplotdata{1,mouse}) ' ' groupkey{group} ' ' columnlabel{action+1}];
                            title(tlt, temptitle, 'fontsize',14);
                            print([graphfolder num2str(tempplotdata{1,mouse}) ' ' groupkey{group} ' ' columnlabel{action+1}], '-dpng');
                        else %will be average col
                            temptitle = [tempplotdata{1,mouse} ' ' groupkey{group} ' ' columnlabel{action+1}];
                            title(tlt,temptitle , 'fontsize',14);
                            print([graphfolder tempplotdata{1,mouse} ' ' groupkey{group} ' ' columnlabel{action+1}], '-dpng');
                        end
                        
                        
                        close all
                        
                    end
                    
                end
                
                
            else %latency plotting
                if row == 5
                    latencytitle = 'ToneToRew';
                    
                else
                    latencytitle = 'RewToRec';
                end
                
                latencybin = latencyedges(2:end);
                
                counts = tempplotdata{row,mouse};
                normcounts = tempplotdata{row+1,mouse};
                
                
                figure
                % tlt=tiledlayout(1,2,'TileSpacing','compact');
                tlt=tiledlayout(1,2);
                nexttile(2)
                imagesc(latencybin, 1, normcounts, [0,1]); % this is the heatmap
                set(gca,'YDir','normal') % put the trial numbers in better order on y-axis
                %fix title
                %title([variables{action} ' Heat Map'],'fontsize',16)
                colormap(gray)
                ylabel('Training Day','fontsize',14)
                xlabel('Time Bin (sec)','fontsize',14)
                cb = colorbar;
                % !Doesn't scale in window properly, fix later
                % ylabel(cb, 'Probability','fontsize',16)
                xticks([0:1:10]);
                
                
                nexttile(1)
                bstack = bar(counts, 'stacked','FaceColor', 'flat');
                % barcolor = colormap(parula(6));
                bincolor = parula(size(counts,2));
                for binnum = 1:size(counts,2)
                    bstack(binnum).CData = bincolor(binnum,:);
                end
                xlabel('Training Day','fontsize',14);
                ylabel('Raw Count','fontsize',14);
                
                
                if ~ischar(tempplotdata{1,mouse})
                    temptitle = [num2str(tempplotdata{1,mouse}) ' ' groupkey{group} ' ' latencytitle];
                    title(tlt, temptitle, 'fontsize',14);
                    print([graphfolder num2str(tempplotdata{1,mouse}) ' ' groupkey{group} ' ' latencytitle], '-dpng');
                    
                else %will be average col
                    temptitle = [tempplotdata{1,mouse} ' ' groupkey{group} ' ' latencytitle];
                    title(tlt, temptitle, 'fontsize',14);
                    print([graphfolder tempplotdata{1,mouse} ' ' groupkey{group} ' ' latencytitle], '-dpng');
                end
                
                
                close all
                
                
            end
            
        end
    end
    
    
end


%% Add group row labels

grouprowlabels = {'Mouse';'Data';'Action Counts'; 'Normalized Action Counts'; 'Tone-Correct Counts'; 'Norm Tone-Correct Counts';'Rew-Rec Counts'; 'Norm Rew-Rec Counts'};

for group = 1:size(groupdata,1)
    
    eval([groupdata{group},'= [grouprowlabels ', groupdata{group} ,' ];']);
    
end


%% Save data in file


xlswrite([outputfolder outputlatfile], AllMeanLat,'AllMeanLat')
xlswrite([outputfolder outputlatfile], AllRewMeanLat,'AllRewMeanLat')
xlswrite([outputfolder outputlatfile], AllAcqTime,'AllAcqTime')

%save all variables together
save([outputfolder outputmatfile]);


%clear data dataname variable row raw

%% Print code version text file

%print the version of the code used
fileID = fopen([graphfolder '\codeused.txt'],'w');
fprintf(fileID, codename);


