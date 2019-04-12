function scores = DPPF_netOV(ActiveDamIndex,nvars,dvars,VarNames)

%NetworkOVcalc is used to calculate the value of ecosystem services that
%rely on interconnected river networks, such as fish habitat and river
%recreation
%ActiveDamIndex=the list of dams that exist for this run and their
%downstream partner(s)
%varargin can contain data on fish habitat and passage, or river
%recreation, or both. No other network-sensitive ecoservices are being
%considered at this time.

%island dam corrections:need to carry over information on the island dams
%to accurately calculate upstream values

if size(ActiveDamIndex,2)>2
    dr=size(ActiveDamIndex,1)-1; %The -1 is for the island data at the end
else
    dr=size(ActiveDamIndex,1);
end

% dud=zeros(dr,1); %debug

if ismember('passage',VarNames)
    if ismember('salmon',VarNames) && (ismember('shad',VarNames) || ismember('alewifeBlueback',VarNames))
        pass=dvars(:,1:2);
        p=2;
    else
        pass=dvars(:,1);
        p=1;
        
    end
else
    p=0;
end

if ismember('riverCanoe',VarNames) && ismember('riverRaft',VarNames)
    r=2;
elseif ismember('riverCanoe',VarNames) || ismember('riverRaft',VarNames)
    r=1;
else
    r=0;
end
netcalc=nvars(1:dr,:);

%Algo to calculate the amount of habitat between each dam, then the FHU
for j=1:dr-1 %the -1 is for the waterhsed at the end
    %   Get the index (indices if island) of downstream dams in order to
    %   subtract value from them
    DownIndex = find(ismember(ActiveDamIndex(1:dr,1),ActiveDamIndex(j,2)));
    if isempty(DownIndex)
        continue
    elseif size(ActiveDamIndex,2)>2 && ActiveDamIndex(j,2)>0%Check for islands. The two most downstream dams on opposing sides of an island overlap in drainage area, therefore need to be corrected by deleting the shared drainage area (and DA-sensitive ecoservices value) from ONE of the island dams.
        if ActiveDamIndex(j,4)>0 %any(netcalc(DownIndex,:)-nvars(j,:)<-0.1) && %&& ActiveDamIndex(j,4)~=IslandCode && sum(ismember(ActiveDamIndex(j,4),ActiveDamIndex(:,4)))>1
            netcalc(j,:)=(nvars(j,:)-nvars(end,:));%new attempt to stop extra fish warning
            netcalc(DownIndex,:)=netcalc(DownIndex,:)-netcalc(j,:);
            %             netcalc(DownIndex,:)=netcalc(DownIndex,:)-(nvars(j,:)-nvars(end,:));
        else
            netcalc(DownIndex,:)=netcalc(DownIndex,:)-nvars(j,:);
        end
    else
        netcalc(DownIndex,:)=netcalc(DownIndex,:)-nvars(j,:);
    end
    
    pt=p;
    rt=r;
    while ActiveDamIndex(DownIndex,1) < 1e9 && pt+rt>0
        if pt
            if sum(pass(j,:)) > 0
                pass(j,:)=pass(j,:).*dvars(DownIndex,1:p);
            else
                pt=0;
            end
        end
        if rt
            dvars(DownIndex,end)=dvars(DownIndex,end)+dvars(j,end-1);
            if r==2
                dvars(DownIndex,end-2)=dvars(DownIndex,end-2)+dvars(j,end-3);
            end
        end
        DownIndex = find(ismember(ActiveDamIndex(1:dr,1),ActiveDamIndex(DownIndex,2)));
        try
        if ActiveDamIndex(DownIndex,1) >= 1e9 && rt
            dvars(DownIndex,end)=dvars(DownIndex,end)+dvars(j,end-1);
            if r==2
                dvars(DownIndex,end-2)=dvars(DownIndex,end-2)+dvars(j,end-3);
            end
            rt=0;
        end
        catch
            dud=1;
        end
    end
end

scores=[];
if any(netcalc<0)
    if any(netcalc<-0.1)
        warning('negative fish warning')
    else
        netcalc(netcalc<0)=0;
    end
end
if p
    if ismember('alewifeBlueback',VarNames) %the alebk dvar is static ie it is not iteratively changed based on downstream/upstream numbers
        netcalc(:,1)=netcalc(:,1) + dvars(1:dr,p+1); %this should be okay to do, because all dam watersheds are now converted to segments between dams (and wshd divides)
    end
    if p==2
        if size(nvars,2)>3
            scores=[scores sum([netcalc(:,1:end-2).*pass(1:dr,1) netcalc(:,end-1).*pass(1:dr,2)],1)];
        else
            scores=[scores sum([netcalc(:,1:end-1).*pass(1:dr,1) netcalc(:,end).*pass(1:dr,2)],1)];
        end
    else
        scores=[scores sum(netcalc(:,1).*pass(1:dr,1))];
    end
end
if r
    %    scores=mean(itemx(itemx(:,end)>0,end),1);
    dvars(dvars(:,end)>1,end)=1;
    netcalc(:,end)=netcalc(:,end).*dvars(1:dr,end);
    scores=[scores max(netcalc(:,end))];
    if r==2
        dvars(dvars(:,end-2)>1,end-2)=1;
        netcalc(:,end-1)=netcalc(:,end-1).*dvars(1:dr,end-2);
        scores=[scores max(netcalc(:,end-1))];
    end
    %     scores=[scores sum(netcalc.*netcalc./sum(netcalc))];
    %     scores=[scores prctile(netcalc(netcalc(:,end)>0,end),90)];
    %     scores=[scores mean(netcalc(netcalc(:,end)>0,end))];
end