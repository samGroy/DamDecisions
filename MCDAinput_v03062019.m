function MCDAinput_v03062019(zipdir,readir,h5file)
tic
%Variables:
%zipdir: the location and name of the zipped folder containing all user
%preference data, produced by Emma's tool
%
%readir: the folder where MOGA-selected scenarios are located
%
%h5file: the h5 file containing all MOGA-selected scenarios and criteria
%values

%if no zipped folder directory given, use default with user preference data
if nargin<1
    zipdir='C:\Users\sgroy\Downloads\Archive.zip';
end
%if no data directory given, use default
if nargin<2
    readir='D:\FoD\PPF\MCDA-PPF\';
    uzdir=[readir 'prefs\'];
end

%get filenames from user preference .csv or .xls files, add to a list
unzip(zipdir,uzdir); %unzip the zipped folder
d1=dir([uzdir '*.csv']);
d2=dir([uzdir '*.xls*']);
d=[d1;d2];
%d=dir([readir '*.xls*']);

%if no h5 file name given, declare h5 file with saved scenarios/data
if nargin<3
    fnm='D:\FoD\PPF\scenarios\';
    hf='MO_12.h5';
    h5file=[fnm hf];
end
x=h5read(h5file,'/Penobscot/x'); %integer returns (decision alternatives) from DamMOGA
f=h5read(h5file,'/Penobscot/f'); %fitness score (criteria quantities) returns from DamMOGA

%save map image directory: maps representing priority decisions will be
%combined with rose plots and copied here
ifnm=[readir 'combo\'];

%declare index list as zeros. Zero entries are omitted prior mapping call.
%Necessary to align results to full dams dataset.
id=zeros(1,length(d));

%declare datavg to take average preference data across all dams and all decision
%alternatives. This was not used in the 2019 test and was commented out.
%datavg=zeros(5,12);

%loop through each user's preference data from AHP sheet, retrieve priority map and rose plot
for i=1:length(d)
    fprintf('reading file # %i\n',i) %report the file number
    try
        dat=xlsread([uzdir d(i).name]); %read xls data
    catch
        fprintf('error reading %s\n',d(i).name) %this occurs if nonumeric entries exist
    end
    if any(any(isnan(dat))) %additional failsafe, if somehow NaN's are in the file.
        fprintf('attempt to import nonumeric values...modifying file #%i\n',i)
        dat(isnan(dat))=0; %%new failsafe, NaNs are very likely to be zeros anyway.
        %continue
    end
    %the following commented out lines were from the original 2018 test. We
    %no longer need to reorder criteria.
    %dat=dat([1 2 7 4 5 3 6],:); %reorder: global, remove, keep, improve turbines, install turbs, improve fish, refurb
    %dat=[dat(1:3,:);mean(dat(4:5,:));dat(6,:);dat(7,:)]; %consolidate: PPF doesn't yet separate improve/install turbines
    %dat(6,:)=mean(dat(4:5,:)); %remove refurb and add double power/fish alternative: PPF doesn't yet incorporate refurb, but allows for overlapping alternatives
    
    %trim user preference data in variable 'dat' down to uniform size,
    %retain all preference data
    dat=dat(1:5,1:12);
    
    %The following variable 'datavg' was used in the original 2018 test, in
    %order to take the average preferences of workshop participants. This
    %functionality was not used in the 2019 test and was taken out.
    %datavg=datavg+dat; %collect all preference data for average after loop complete
    
    %call external function prefDA:
    %use fn prefDA to generate weights based on global weights and local
    %(decision alternative) weights based on frequency of decision
    %alternatives in each scenario. Produces a unique set of preference
    %weights for each decision scenario (consisting of a set of dam
    %decisions). This essentially collapses the different preferences of
    %each decision alternative, providing a weighted average preference
    %value for each criteria.
    w=prefDA(dat,x); 
    
    %call external function MultiRank
    %Use MultiRank to generate ranked list of scenarios, given the weighted
    %preferences of the user (and modified in prefDA to collapse across
    %decision alternatives)
    rnk=MultiRank(f,w,'normax','weightedSum'); %rank scenarios by weights. 'normax' means normalize relative to maximum value. 'weightedSum' means use the weighted sum approach
    
    %In 2018 test, I created rose plots in real time. For the 2019 test,
    %all possible scenarios are collected and rose plots are pre-generated,
    %so this has been commented out.
    %RosePlots(f,1,w,[0 0 0],[ifnm d(i).name(1:2)],'weightedProduct','normax')
    %close(gcf)
    
    %The highest rank scenario is the 'best' scenario. Select it.
    id(i)=rnk(1);
    
    %retrieve the map that corresponds to the highest ranked scenario
    mapdir=dir([readir 'combo\*_' sprintf('%i.png',id(i)-1)]);
    %Rename the map based on the user filename and scenario number.
    %copy the map to the save folder
    cpyname=[extractBefore(d(i).name,'.') sprintf('__%i.png',id(i)-1)];
    copyfile([readir '\combo\' mapdir.name()],[readir '\selected\' cpyname]);
end

%The following commented out lines were used for 1) averaging preference data
%and 2) creating scenario maps from scratch, these functions were not used in the 2019 test.
% %calculate average weights and select top 3 scenarios based on averages
% datavg=datavg./length(d);
% w=prefDA(datavg,x); %use fn prefDA to generate weights based on global weights, local weights, and frequency of decision alternatives
% rnk=MultiRank(f,w,'normax','weightedProduct'); %rank scenarios by weights
% id=[id rnk(1:3)'];
% RosePlots(f,1,w,[0 0 0],[ifnm 'AvgRank1'],'weightedProduct','normax')
% close(gcf)
% RosePlots(f,2,w,[0 0 0],[ifnm 'AvgRank2'],'weightedProduct','normax')
% close(gcf)
% RosePlots(f,3,w,[0 0 0],[ifnm 'AvgRank3'],'weightedProduct','normax')
% close(gcf)
% %%quick debug only
% %id=unique(id);
% 
% id=id-1;
% 
% cd(fnm)
% system(['ScenarioDamSelector04.py ' fnm hf ' Penobscot MO 07 ' sprintf('%i,',id)]);
toc %check solution speed