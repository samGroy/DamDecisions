function scores = Barr_netOV(ActiveDamIndex,nvars,pv,av,rv)

%NetworkOVcalc is used to calculate the value of ecosystem services that
%rely on interconnected river networks, such as fish habitat and river
%recreation
%ActiveDamIndex=the list of dams that exist for this run and their
%downstream partner(s)
%varargin can contain data on fish habitat and pvage, or river
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

% if ismember('pvage',VarNames)
%     if ismember('salmon',VarNames) && (ismember('shad',VarNames) || ismember('alewifeBlueback',VarNames))
%         pv=dvars(:,1:2);
%         p=2;
%     else
%         pv=dvars(:,1);
%         p=1;
%
%     end
% else
%     p=0;
% end
%
% if ismember('riverCanoe',VarNames) && ismember('riverRaft',VarNames)
%     r=2;
% elseif ismember('riverCanoe',VarNames) || ismember('riverRaft',VarNames)
%     r=1;
% else
%     r=0;
% end
if ~isnan(pv(1,1))
    p=1;
else
    p=0;
end
if size(nvars,2)>size(pv,2)
    r=1;
else
    r=0;
end
netcalc=nvars(1:dr,:);
pass=pv;
%Algo to calculate the amount of habitat between each dam, then the FHU
for j=1:dr-1 %the -1 is for the waterhsed at the end *******WATCH THAT, there can be more than 1 watershed!!
    %   Get the index (indices if island) of downstream dams in order to
    %   subtract value from them
    DownIndex = find(ismember(ActiveDamIndex(1:dr,1),ActiveDamIndex(j,2)));
    if isempty(DownIndex)
        continue
    elseif size(ActiveDamIndex,2)>2 && ActiveDamIndex(j,2)>0%Check for islands. The two most downstream dams on opposing sides of an island overlap in drainage area, therefore need to be corrected by deleting the shared drainage area (and DA-sensitive ecoservices value) from ONE of the island dams.
        if ActiveDamIndex(j,4)>0 %any(netcalc(DownIndex,:)-nv(j,:)<-0.1) && %&& ActiveDamIndex(j,4)~=IslandCode && sum(ismember(ActiveDamIndex(j,4),ActiveDamIndex(:,4)))>1
            netcalc(j,:)=(nvars(j,:)-nvars(end,:));%new attempt to stop extra fish warning
            netcalc(DownIndex,:)=netcalc(DownIndex,:)-netcalc(j,:);
            %             netcalc(DownIndex,:)=netcalc(DownIndex,:)-(nv(j,:)-nv(end,:));
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
            if sum(pv(j,:)) > 0
                pv(j,:)=pv(j,:).*pass(DownIndex,:);
            else
                pt=0;
            end
        end
        if rt
            rv(DownIndex,end)=rv(DownIndex,end)+rv(j,end-1);
        end
        DownIndex = find(ismember(ActiveDamIndex(1:dr,1),ActiveDamIndex(DownIndex,2)));
        try
            if ActiveDamIndex(DownIndex,1) >= 1e9 && rt
                rv(DownIndex,end)=rv(DownIndex,end)+rv(j,end-1);
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
    if ~isnan(av) %the alebk dvar is static ie it is not iteratively changed based on downstream/upstream numbers
        netcalc(:,1)=netcalc(:,1) + av(1:dr); %this should be okay to do, because all dam watersheds are now converted to segments between dams (and wshd divides)
    end
    scores=[scores sum(netcalc(:,1:size(pv,2)).*pv(1:dr,:),1)]; %product of MFHUs and pv's
end
if r
    %    scores=mean(itemx(itemx(:,end)>0,end),1);
    rv(rv(:,end)>1,end)=1;
    netcalc(:,end)=netcalc(:,end).*rv(1:dr,end);
    scores=[scores max(netcalc(:,end))];
end