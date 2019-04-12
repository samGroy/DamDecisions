function [BarrIndex, v] = ...
    BarrOpt_prep(watershed, basedir, varargin) %%%%%,costRemove, lakeRecreationRm, alewifeBluebackRm

%9-5-18
%Code built from DPPF suite to incorporate culvert metrics and lump old
%metrics
%watershed = name of watershed, must be in your BarrIndex database
%basedir = base directory for accessing barrier index and data
%varargin = additional arguments used to refine/filter the returned database

%watershed = 'Connecticut';
BarrDataF = sprintf('%s/BarrData/',basedir);
BarrIndexF = sprintf('%s/BarrIndex/',basedir);

D=readtable(sprintf('%sBarr_4-9-19.csv',BarrDataF));

%%%
%Order of data:
%identifiers (BID, IsBarrier, BarTyp)
%infrastructure failure risk (OT_RISK)
%diadromous fish (BkTHabA, BKT_P, RhHabA, AlRm5d_skm, Al5Dam, AleP, BbkP)
%dam utilities (Cpcty_kW, RsvrVolEst)
%recreational boating (FRRU, fdy300, fdr300cfs, RsvrRm_skm**replace with water quality corrected form**)
%decision cost (DnCost)
%%%

%criteria: names of the outputs, criteria used for analysis
criteria={'Infrastructure Failure Risk',...
    'Diadromous Fish Habitat',...
    'Dam Utilities',...
    'Recreational Boating',...
    'Decision Costs'};

%fields
id_fields=cell(1,1);
data_fields=cell(1,5);
%IDs
id_fields{1}={'BID','IsBarrier','BarTyp','HUC8','HUC12','GEOCODE'};

%crop data table to the data_fields for this analysis, considering the
%variables of interest provided by user
%also crop the criteria
if any(strcmp(inputname(nargin),{'combo','variables','vars','ovars','ov','ovs','data_fields'}))
    if iscell(varargin{nargin-2})
        criteria_codes=[];
        %Infrastructure Failure Risk
        if any(contains(varargin{nargin-2},'risk','IgnoreCase',true) ...
                | contains(varargin{nargin-2},'infrastructure','IgnoreCase',true) ...
                | strcmp(varargin{nargin-2},'I'))
            
            criteria_codes=[criteria_codes 1];
            data_fields{1}=[data_fields{1} {'OT_RISK'}];
        end
        %lots of fish names
        if any(contains(varargin{nargin-2},'fish','IgnoreCase',true) ...
                | contains(varargin{nargin-2},'diadromous','IgnoreCase',true) ...
                | contains(varargin{nargin-2},'habitat','IgnoreCase',true) ...
                | strcmp(varargin{nargin-2},'F'))
            
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'RhHabA','RH_P','AlRm_5d_sm','BkTHabA','BKT_P',...
                'SalmoHabA','ASAL_P','ShadHabA','ASHA_P','SLamHabA','SL_P','EelHabA','AE_P'}];
        end
        if any(contains(varargin{nargin-2},'raw habitat','IgnoreCase',true))
            
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'RhHabA','RH_P','BkTHabA','BKT_P',...
                'SalmoHabA','ASAL_P','ShadHabA','ASHA_P','SLamHabA','SL_P','EelHabA','AE_P'}];
            
            %If user wants to do a "#habitat vs #barriers" comparison
            %rather than use passage probs
            if any(contains(varargin{nargin-2},'number of barriers','IgnoreCase',true))
                D.RH_P=floor(D.RH_P);
                D.BKT_P=D.RH_P;
                D.ASAL_P=D.RH_P;
                D.ASHA_P=D.RH_P;
                D.SL_P=D.RH_P;
                D.AE_P=D.RH_P;
            end
        end
        %specify river herring
        if any(contains(varargin{nargin-2},'herring','IgnoreCase',true)) ...
                || any(contains(varargin{nargin-2},'alewife','IgnoreCase',true))
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'RhHabA','AlRm_5d_sm','RH_P'}];
        end
        %specify trout
        if any(contains(varargin{nargin-2},'trout','IgnoreCase',true))
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'BkTHabA','BKT_P'}];
        end
        %specify salmon
        if any(contains(varargin{nargin-2},'salmon','IgnoreCase',true))
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'SalmoHabA','ASAL_P'}];
        end
        %specify shad
        if any(contains(varargin{nargin-2},'shad','IgnoreCase',true))
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'ShadHabA','ASHA_P'}];
        end
        %specify eel
        if any(contains(varargin{nargin-2},'eel','IgnoreCase',true))
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'EelHabA','AE_P'}];
        end
        %specify sea lamprey
        if any(contains(varargin{nargin-2},'lamprey','IgnoreCase',true))
            criteria_codes=[criteria_codes 2];
            data_fields{2}=[data_fields{2} {'SLamHabA','SL_P'}];
        end
        %Dam Utilities
        if any(contains(varargin{nargin-2},'dam','IgnoreCase',true) ...
                | contains(varargin{nargin-2},'utilities','IgnoreCase',true) ...
                | strcmp(varargin{nargin-2},'U'))
            
            criteria_codes=[criteria_codes 3];
            data_fields{3}=[data_fields{3} {'Cpcty_kW','RsvrVolEst'}];
        end
        %specify power
        if any(contains(varargin{nargin-2},'power','IgnoreCase',true))
            criteria_codes=[criteria_codes 3];
            data_fields{3}=[data_fields{3} {'Cpcty_kW'}];
        end
        %specify storage
        if any(contains(varargin{nargin-2},'storage','IgnoreCase',true))
            criteria_codes=[criteria_codes 3];
            data_fields{3}=[data_fields{3} {'RsvrVolEst'}];
        end
        %Recreational Boating
        if any(contains(varargin{nargin-2},'rec','IgnoreCase',true) ...
                | strcmp(varargin{nargin-2},'R'))
            
            criteria_codes=[criteria_codes 4];
            data_fields{4}=[data_fields{4} {'FRRU','fdr300cfs','fdy300','FLRU'}];
        end
        if any(contains(varargin{nargin-2},'river','IgnoreCase',true))
            criteria_codes=[criteria_codes 4];
            data_fields{4}=[data_fields{4} {'FRRU','fdr300cfs','fdy300'}];
        end
        if any(contains(varargin{nargin-2},'lake','IgnoreCase',true))
            criteria_codes=[criteria_codes 4];
            data_fields{4}=[data_fields{4} {'FLRU'}];
        end
        
        %Decision Cost
        if any(contains(varargin{nargin-2},'cost','IgnoreCase',true) ...
                | strcmp(varargin{nargin-2},'C'))
            
            criteria_codes=[criteria_codes 5];
            data_fields{5}=[data_fields{5} {'DnCost'}];
        end
        criteria_codes=unique(floor(criteria_codes));
        
    else
        %Infrastructure Failure Risk
        data_fields{1}={'OT_RISK'};
        %Diadromous Fish
        data_fields{2}={'RhHabA','AlRm_5d_sm','RH_P','BkTHabA','BKT_P'};
        %Dam Utilities
        data_fields{3}={'Cpcty_kW','RsvrVolEst'};
        %Recreational Boating
        data_fields{4}={'FRRU','fdr300cfs','fdy300','FLRU'};%'RsvrRm_skm',
        %Decision Cost
        data_fields{5}={'DnCost'};
        
        criteria_codes=unique(floor(varargin{nargin-2}));
    end
    data_fields=data_fields(criteria_codes);
    criteria=criteria(criteria_codes);
end
D=D(:,[id_fields{:} data_fields{:}]);

%if all watersheds in New England are used for assessment, append them in
%one dataset, append DamIndex
if iscell(watershed)
    BarrIndex=[];
    for i=1:length(watershed)
        di=dlmread(sprintf('%s%sDownDamList.txt',BarrIndexF,watershed{i}),',');
        if size(di,2)<4
            di(:,3:4)=0;
        end
        BarrIndex=[BarrIndex; di];
    end
    if sum(sum(BarrIndex(:,3:4)))==0
        BarrIndex=BarrIndex(:,1:2);
    end
elseif strcmpi(watershed,'All') || ~isnan(str2double(watershed(7:end))) %the second check is if a scale filter is applied that would require multiple watersheds
    d=dir([BarrIndexF '*' 'DownDamList.txt']);
    BarrIndex=[];
    for i=1:length(d)
        di=dlmread([d(i).folder '/' d(i).name],',');
        if strcmp(d(i).name,'PenobscotLowerDownDamList.txt')
            continue
        end
        if size(di,2)<4
            di(:,3:4)=0;
        end
        BarrIndex=[BarrIndex; di];
    end
else
    BarrIndex=dlmread(sprintf('%s%sDownDamList.txt',BarrIndexF,watershed),',');
end

%First, filter to keep only dams in the watershed of interest
[~,ind1]=ismember(D{:,1},BarrIndex(:,1)); %BID is always first column in table
D=D(ind1>0,:);
ind2=ind1(ind1>0,:);
BarrIndex=BarrIndex(ind2,:);
%DamDat=DamDat(ind2,:);

%find the index of the watershed, used later to make sure watershed ends up
%on the bottom of the list
% indW=find(DamIndex(:,1)>=1e9);
Wid=BarrIndex(BarrIndex(:,1)>=1e9,:);
Wdat=D(BarrIndex(:,1)>=1e9,:);

%if islands exist with dams on both sides
if size(BarrIndex,2)>2
    Iid=unique(BarrIndex(BarrIndex(:,4)>0,4)); %hardcoded for Milford, the island hub dam
    Idat=D(BarrIndex(:,1)==Iid,:);
end

%sub-selection based on traits, e.g., hydropower dams, reservoir dams, etc.
%First, did user request removal of 'zero' barriers (those that don't
%provide or impact any measurable value.)
for i=1:length(data_fields)
    if contains(data_fields{i},'DnCost')
        continue
    end
    select_fields{i}=data_fields{i}(~contains(data_fields{i},'_P'));
end
if any(strcmp('zero',varargin))
    x=max(table2array(D(:,[select_fields{:}])),[],2);
    x(x>0)=1;
    x(x<0)=0;
elseif any(strcmp('percentile',varargin))
    if isnumeric(varargin{find(strcmp(varargin,'percentile'))+1})
        if size(varargin{find(strcmp(varargin,'percentile'))+1},2)==1
            percnum=varargin{find(strcmp(varargin,'percentile'))+1};
        else
            percnum=10;
        end
    end
    p_threshold=prctile(table2array(D(:,[select_fields{:}])),percnum,1);
    x=sum(table2array(D(:,[select_fields{:}]))-p_threshold,2);
    x(x>=0)=1;
    x(x<0)=0;
else
    x=ones(length(BarrIndex),1);
end

%default remove anything that isn't considered a "barrier"
x(D{:,'IsBarrier'}==0)=0;

%sub-sub-selection based on scale. If user fed scalevar to this function, use
%the input to select only barreris within the region defined by the
%scalevar code. scalevar = 'scale-code', for example, 'HUC08-1060001', 'HUC12-10600010306
%or 'TOWNS-5070'
if ~isnan(str2double(watershed(7:end)))
    if strcmp(watershed(1:6),'HUC08-')
        x(D{:,'HUC8'}~=str2double(watershed(7:end)))=0;
    elseif strcmp(watershed(1:6),'HUC12-')
        x(D{:,'HUC12'}~=str2double(watershed(7:end)))=0;
    elseif strcmp(watershed(1:6),'TOWNS-')
        x(D{:,'GEOCODE'}~=str2double(watershed(7:end)))=0;
    end
end

%user specify just culverts or dams?
if any(strcmp('culvs',varargin) | strcmp('culverts',varargin))
    x(D{:,'BarTyp'}~=2)=0;%culverts only
elseif any(strcmp('dams',varargin)) %dams only
    x(D{:,'BarTyp'}~=1)=0;
end

%consider just current infrastructure e.g. no new dams?
if any(strcmp('crh1',varargin))
    x(D{:,'BarTyp'}==3)=0;
end
%Make sure watershed shows up at the end
D(BarrIndex(:,1)>=1e9,:)=[];
x(BarrIndex(:,1)>=1e9,:)=[];
BarrIndex(BarrIndex(:,1)>=1e9,:)=[];
BarrIndex=[BarrIndex; Wid];
D=[D; Wdat];
x=[x; ones(size(Wid,1),1)];

%%Output variable definitions
%trim the index and data lists
NewIndex=DPPF_idx(x,BarrIndex,0);
D=D(x==1,:);
BarrIndex=NewIndex;

%if there's an island and the island hub was removed during selection
if size(BarrIndex,2)>2
    if sum(sum(BarrIndex(:,3:4)))>0
        BarrIndex(end+1,1)=Iid;
        D(end+1,:)=Idat;
        x(end+1)=1;
    else
        BarrIndex=BarrIndex(:,1:2);
    end
end

%organize output variable with data and criterion names
v=cell(2,length(criteria));
for i=1:length(criteria)
    v{1,i}=criteria{i};
    v{2,i}=D(:,data_fields{i});
end
