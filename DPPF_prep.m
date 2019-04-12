function [DamIndex, v] = ...
    DPPF_prep(watershed, basedir, varargin) %%%%%,costRemove, lakeRecreationRm, alewifeBluebackRm
%DamIndexSelector
%Use this script to select the watershed/dams you want to PPF
%WARNING: this function is messy as hell

% %make sure args are in char format
% if any(strcmp('percentile',varargin)) && ~any(strcmp('batch',varargin))
%     if isnumeric(varargin{find(strcmp(varargin,'percentile'))+1})
%         percnum=varargin{find(strcmp(varargin,'percentile'))+1};
%         varargin{find(strcmp(varargin,'percentile'))+1}=num2str(varargin{find(strcmp(varargin,'percentile'))+1});
%     else
%         percnum=10;
%     end
% end

%varargin: list of filters for dam selection
%watershed = 'Connecticut';
DamDataF = sprintf('%s/DamData/',basedir);
DamIndexF = sprintf('%s/DamIndex/',basedir);

if strcmp(watershed,'PenobExpCap')
    DamDat=dlmread(sprintf('%sNEdamsOut8-5-17clipExpCap.txt',DamDataF),'\t',1,0);
    PowerExp=34;
    watershed='Penobscot';
elseif strcmp(watershed,'PenobFishExp')
    DamDat=dlmread(sprintf('%sNEdamsOut8-5-17clipExpCap.txt',DamDataF),'\t',1,0);
    FishPassExp=27;
    watershed='Penobscot';
elseif strcmp(watershed,'PenobFishExpCap')
    DamDat=dlmread(sprintf('%sNEdamsOut8-5-17clipExpCap.txt',DamDataF),'\t',1,0);
    PowerExp=34;
    FishPassExp=27;
    watershed='Penobscot';
else
%     D=readtable(sprintf('%sNEdamsOut6-12-18.txt',DamDataF));
    D=readtable(sprintf('%sNEdamsOut2-5-19.txt',DamDataF));
    DamDat=[D.Cpcty_kW ...%1
        D.FERC_ID ...%2
        D.DID ...%3
        D.CRH ...%4
        D.Resvr_skm ...%5
        D.Hazard ...%6
        D.PassUp ...%7
        D.PassDn ...%8
        D.PassUpSal ...%9
        D.PassDnSal ...%10
        D.DA_skm ...%11
        D.FRRU ...%12 Replaces FastWtr_km
        D.StrcRD200 ...%13
        D.StrcRD1000 ...%14
        D.PopServ ...%15
        D.PWSID ...%16
        D.Omit ...%17
        D.SUM_AtlSal ...%18
        D.SUM_ShadMF ...%19
        D.SUM_SalNOA ...%20
        D.Q0001E ...%21
        D.SUM_AlBkMF ...%22 %RH
        D.RsvrVolEst ...%23
        D.fdr150cfs ...%24 %D.ResQrat ...%24... %RiverRec
        D.fdr300cfs ...%25
        D.fdr2000cfs ...%26
        D.ShLnDamAdd...D.RsvrRm_skm ... %27 (25) LakeRec
        D.N_rm_DLoss ... %28 (26) Nitro
        D.AlRm5d_skm ...%29 (27) RH
        D.Al5Dam ...%30 (28) RH 
        D.fdy150cfs ... %31 %D.FracSeason ... %(29) RiverRec
        D.fdy300cfs ... %32
        D.fdy2000cfs ... %33
        D.ExpPow ... %34 (30) %power alternatives
        D.TotalCostE ...%35 (31) %cost
        D.ExpCost];%36 (32) %alternatives
    if strcmp(watershed,'PenobscotLower')
        DamDat(:,18)=D.SUM_AtlS_1;
%         DamDat(:,22)=D.SUM_AlBk_1;
%         DamDat(:,30)=D.Al5DamLP;
        DamDat(DamDat(:,3)==(1e9+1),27)=D.SLDA_LowPe(D.SLDA_LowPe>0);
    end
    D=[];
%     DamDat=dlmread(sprintf('%sNEdamsOut10-13-17clip.txt',DamDataF),'\t',1,0);
%     DamDat=dlmread('\\NEdamsOut7-19-17clip.txt','\t',1,0);
end

%if all watersheds in New England are used for assessment, append them in
%one dataset, append DamIndex
if strcmp(watershed,'all')
    d=dir([DamIndexF '*' 'DownDamList.txt']);
    DamIndex=[];
    for i=1:length(d)
        di=dlmread([d(i).folder '/' d(i).name],',');
        if strcmp(d(i).name,'PenobscotLowerDownDamList.txt')% || strcmp(d(i).name,'CoastalDownDamList.txt')
            continue
        end
        if size(di,2)<4
            di(:,3:4)=0;
        end
        DamIndex=[DamIndex; di];
    end
else
    DamIndex=dlmread(sprintf('%s%sDownDamList.txt',DamIndexF,watershed),',');
end

%filter to keep only dams in the watershed of interest
[~,ind1]=ismember(DamDat(:,3),DamIndex(:,1)); 
DamDat=DamDat(ind1>0,:);
ind2=ind1(ind1>0,:);
DamIndex=DamIndex(ind2,:);
%DamDat=DamDat(ind2,:);

%dam dependencies: drinking water dams. Need to make sure what dams in a
%network will become absent
drinktrig=0;
if length(DamDat(DamDat(:,16)>0,16))~=length(unique(DamDat(DamDat(:,16)>0,16)))
    DrinkIndex=find(DamDat(:,16)>0);
    DamDat(:,end+1)=0;
    drinktrig=1;
    for i=1:length(DrinkIndex)
        DrinkDamCount=sum(DamDat(:,16)==DamDat(DrinkIndex(i),16));
        DamDat(DrinkIndex(i),end)=DrinkDamCount;
    end
end
%make sure DamIndex values line up before selecting the dams you want. Keep
%track of the index for the watershed.
%ind2=find(DamIndex(:,1)==DamDat(:,3));
%DamIndex=DamIndex(ind2,:);

%find the index of the watershed, used later to make sure watershed ends up
%on the bottom of the list
% indW=find(DamIndex(:,1)>=1e9);
Wid=DamIndex(DamIndex(:,1)>=1e9,:);
Wdat=DamDat(DamIndex(:,1)>=1e9,:);

%if islands exist with dams on both sides
if size(DamIndex,2)>2
    Iid=unique(DamIndex(DamIndex(:,4)>0,4)); %hardcoded for Milford, the island hub dam
    Idat=DamDat(DamIndex(:,1)==Iid,:);
end

%sub-selection based on traits, e.g., hydropower dams, reservoir dams, etc.
%ArgCount = nargin - 1; %one argument to select the watershed, the rest are to filter dam selection
%if ArgCount
x=ones(size(DamIndex,1),1); %REMEMBER this x is totally independent of the x that will be used in NSGA-II!!!!!
idx=0;

%did user provide codes to services of interest
for i=1:length(varargin)
    if isnumeric(varargin{i})
        %             if length(varargin{i})>1
        idx=i;
        %             end
    end
    if i==length(varargin) && ~idx
        varargin{end+1}=1:1:13;
    end
end

if ismember(1,varargin{end})
    x(DamDat(:,1)>0)=x(DamDat(:,1)>0)+1;%hydropower
end
if ismember(2,varargin{end})
    x(DamDat(:,15)>0)=x(DamDat(:,15)>0)+1;%drink
end
if ismember(3,varargin{end})
    x(DamDat(:,23)>0)=x(DamDat(:,23)>0)+1;%storage
end
if ismember(4,varargin{end})
    x(DamDat(:,28)>0)=x(DamDat(:,28)>0)+1;%nitrogen
end
if ismember(5,varargin{end})
    x(DamDat(:,13)>0)=x(DamDat(:,13)>0)+1;%property
end
if ismember(13,varargin{end})
    x(DamDat(:,27)>0)=x(DamDat(:,27)>0)+1;%lakeRecreation
end
if ismember(8,varargin{end})
    x(sum(DamDat(:,[22 30 29]),2)>0)=x(sum(DamDat(:,[22 30 29]),2)>0)+1;%alewifeBlueback
end
if ismember(9,varargin{end})
    x(DamDat(:,19)>0)=x(DamDat(:,19)>0)+1;%shad
end
if ismember(10,varargin{end})
    x(DamDat(:,18)>0)=x(DamDat(:,18)>0)+1;%salmon
end
if ismember(12,varargin{end})
    x(sum(DamDat(:,[12 25 32]),2)>0)=x(sum(DamDat(:,[12 25 32]),2)>0)+1;%riverRec
end
if ismember(7,varargin{end})
    x(DamDat(:,6)>1)=x(DamDat(:,6)>1)+1;%hazard
end
if ismember(6,varargin{end})
    x(DamDat(:,35)>0)=x(DamDat(:,35)>0)+1;%cost of removal
end
if any(strcmp('crh1',varargin))
    x(DamDat(:,4)~=1)=0;%present dams only
elseif any(strcmp('crh2',varargin))
    x(~ismember(DamDat(:,4),[1,2]))=0;
end


%Binarize the dams that were selected
if any(strcmp('zero',varargin))
    x(x==1)=0;
    x(x>1)=1;
    %remove omits: the dams that have too small a watershed
    x(DamDat(:,17)>0)=0;
elseif any(strcmp('percentile',varargin))
    if isnumeric(varargin{find(strcmp(varargin,'percentile'))+1})
        if size(varargin{find(strcmp(varargin,'percentile'))+1},2)==1
            percnum=varargin{find(strcmp(varargin,'percentile'))+1};
        else
            percnum=10;
        end
    end
    didx=[1 15 23 27 28 13 30 22 19 18 12 6 35]; %maps out where each v column is located in DamData...
    dd=DamDat(:,didx(varargin{idx}(1:end-1)));
%     ddz=prctile(dd(x>0,:),percnum,1);
    for i=1:size(dd,2)
        ddz(i)=prctile(dd(x>0 & dd(:,i)>0,i),percnum);
    end
    for i=1:size(dd,2)
        dd(dd(:,i)<=ddz(i),i)=0;
        dd(dd(:,i)>ddz(i),i)=1;
    end
    dd=sum(dd,2);
    thresh=0;
    x(dd<=thresh)=0;
    x(dd>thresh & x>0)=1;
    %remove omits: the dams that have too small a watershed
    x(DamDat(:,17)>0)=0;
else
    x(x>1)=1;
%     x(x<max(x))=0;
%     x(x==max(x))=1;
end

%Make sure watershed shows up at the end
DamDat(DamIndex(:,1)>=1e9,:)=[];
x(DamIndex(:,1)>=1e9,:)=[];
DamIndex(DamIndex(:,1)>=1e9,:)=[];
DamIndex=[DamIndex; Wid];
DamDat=[DamDat; Wdat];
x=[x; ones(size(Wid,1),1)];

%if there's an island and the island hub was removed during selection
if size(DamIndex,2)>2
    DamIndex(end+1,1)=Iid;
    DamDat(end+1,:)=Idat;
    x(end+1)=1;
end

%trim the lists
ActiveDamIndex=DPPF_idx(x,DamIndex,0);
DamDat=DamDat(x==1,:);
DamIndex=ActiveDamIndex;
%end
%write out the different variables used for indexing and OV calc
%indexy stuff
% DID=DamDat(:,3);
% crh=DamDat(:,4);
% DA=DamDat(:,11);
% omit=DamDat(:,17);
v=cell(2,13);

%hydropower
if any(strcmp('ExpPow',varargin))
    power=[DamDat(:,1) DamDat(:,34)];
    power(DamDat(:,4)==2,:)=power(DamDat(:,4)==2,[2 1]);
else
    power=DamDat(:,1);
end
v{1,1}='power';
v{2,1}=power;
powerDependency=DamDat(:,2);
% power=[kw powerDependency];
%power=kw; %for now, forget dependencies besides the major one in Milford

%drinking water
drink=DamDat(:,15);
if drinktrig
    drink=[drink DamDat(:,16) DamDat(:,end)];
end
v{1,2}='drink';
v{2,2}=drink;

%storage
storage=DamDat(:,23);
v{1,3}='storage';
v{2,3}=storage;

%nitrogen
nitrogen=DamDat(:,28);
v{1,4}='nitrogen';
v{2,4}=nitrogen;

%landowner stuff
Struc200=DamDat(:,13);
Struc1000=DamDat(:,14);
property=Struc200;%[Struc200 Struc1000];
v{1,5}='property';
v{2,5}=property;

%cost
%dam heights are DamDat(:,26)
%costRemove=DamDat(:,27); %canceled for now, but get it sorted
cost=DamDat(:,35);
if any(strcmp('ExpPow',varargin))
    cost=[cost DamDat(:,36)];
end
if any(strcmp('ExpFish',varargin))
    cost=[cost cost(:,1).*2]; %simple rule: fish passage is twice as expensive as dam removal
end
if any(strcmp('ExpPow',varargin))
    cost(DamDat(:,4)==2,1)=0; %if ExpPow, crh2 dams counted as expansion dams, so there is no removal fee, only construction fee. If not ExpPow, assume these dams already exist.
end
v{1,6}='cost';
v{2,6}=cost;

%hazard
hazard=DamDat(:,6);
hazard(hazard>0)=hazard(hazard>0)-1; %assumes that if state hasn't rated hazard for the dam (0), it is probably low risk. Assumes they use some prioritization scheme to target high hazard dams first.
%hazard(hazard<3)=0; 
%hazard(hazard>0)=1;%fitness function trick to sum the highest hazard dams. This minimizes the calculations needed to do that, by removing lower hazard dams.
v{1,7}='hazard';
v{2,7}=hazard;

%fish
%convert all to km^2 from m^2, alewives are already in those units
alewifeBlueback=DamDat(:,30)+DamDat(:,22)./1e6;
alewifeBlueback=[alewifeBlueback DamDat(:,29)];
v{1,8}='alewifeBlueback';
v{2,8}=alewifeBlueback;

% ales=DamDat(:,26);
shad=DamDat(:,19)./1e6;
v{1,9}='shad';
v{2,9}=shad;

pass=DamDat(:,7).*DamDat(:,8);
%alewifeBlueback=DamDat(:,22)./1e6+alewife;
%alewifeBluebackRm=DamDat(:,25)./1e6+DamDat(:,27);

salmon=DamDat(:,18)./1e6;
v{1,10}='salmon';
v{2,10}=salmon;

salmonNOAA=DamDat(:,20)./1e6;
Spass=DamDat(:,9).*DamDat(:,10);
%fish=[Alewife Blueback Shad SalmonNOAA]; %let user build this manually.
%Not perfect, but allows for fewer OVs if desired.
passage=[pass Spass];
if any(strcmp('ExpFish',varargin))
    passage=[passage pass+((1-pass)*.5) Spass+((1-Spass)*.5)];
end
v{1,11}='passage';
v{2,11}=passage;

%River recreation by canoe, convert FRRU m^2 to km^2
riverRecreation=[DamDat(:,12)./1e6 DamDat(:,25) DamDat(:,32)];
% riverRecreation(riverRecreation(:,3)==0,3)=0.333;
v{1,12}='riverCanoe';
v{2,12}=riverRecreation;

%River recreation by canoe, convert FRRU m^2 to km^2
riverRecreation=[DamDat(:,12)./1e6 DamDat(:,26) DamDat(:,33)];
% riverRecreation(riverRecreation(:,3)==0,3)=0.333;
v{1,13}='riverRaft';
v{2,13}=riverRecreation;

% Aggregate recreation
%recreation = [lakeRecreation riverRecreation];

%flatwater rec
lakeRecreation=DamDat(:,27);
v{1,14}='lakeRecreation';
v{2,14}=lakeRecreation;

%discharge: model mean discharge in m^3s^-1
discharge=DamDat(:,21);

%ag: constants that are used to standardize and aggregate similar ES
%ag1-3: carrying capacity of river herring, shad, salmon
%ag4: drainage area of the watershed in consideration
%ag5-6: New England lake and river rec densities, respectively
ag=[21.3 24.216 3.002 0.0335 0.0275 sum(DamDat(DamDat(:,3)>=1e9,11)) 2 2]; %ag(7) triggers aggregation (1) or not (0)
if any(strcmp('ExpPow',varargin))
    ag(8)=ag(8)+1;
end
if any(strcmp('ExpFish',varargin))
    ag(8)=ag(8)+1;
end
if any(strcmp('ExpPow',varargin)) && any(strcmp('ExpFish',varargin))
    ag(8)=ag(8)+1;
end
v{1,15}='ag';
v{2,15}=ag;
%ag(7) = 2 default sum canoe and rafting FRRU

%PowerExpansion quoted on Penobscot
if exist('PowerExp','var')
    power=power+DamDat(:,PowerExp);
end
if exist('FishPassExp','var')
    passage=[passage(:,1)+DamDat(:,FishPassExp) passage(:,2)+DamDat(:,FishPassExp+1)];
    passage(passage>0 & passage<0.75)=0.75;
end
