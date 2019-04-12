function ActiveDamIndex=DPPF_idx(x,DamIndex,calltype)

%DownDamReorganizer: use to remove dams from the downstream dam list index
%that are removed by NSGA-II, and therefore not required for the objective
%variable calculations invoked in DPPF_fitfn.m

%DamIndex = the original, static downstream dam index defined during
%preprocessing

%x = the list of dams that are kept or removed, produced by NSGA-II
%The watershed ID is pinned to the end but is not picked by NSGA-II
%If there are islands, the watershed ID will be followed by those data. The
%island values are not included in x.

%ActiveDamList = the downstream dam list index defined for this cycle,
%updated from the NSGA-II removals

%array the removed dams. Island data on the end (if it exists) is not
%included in this
rmDams = DamIndex(x==0,1);

%flag for loop while searching for downstream dams that still exist
gone = 1;
if size(DamIndex,2)>2
    island=1;
    cols=2:3;
    %     ActiveDamIndex = DamIndex(1:end-1,:); %assume there's only one island, therefore only one extra tag after the watershed tag
else
    island=0;
    cols=2;
    %     ActiveDamIndex = DamIndex(1:end,:); %just the watershed tag is removed to sync with the x array
end
ActiveDamIndex=DamIndex;
prevDownDam=0;
%loop to update active dam list
for k=cols %for islands
    for i=1:length(rmDams)
        updateDams = DamIndex(DamIndex(1:end-island,k)==rmDams(i),1); %grab dams w/ removed downstream dams to update their new downstream neighbor
        updateDams = updateDams(~ismember(updateDams,rmDams)==1); %some of these dams will themselves be removed, ignore those to be removed from the list at the end
        if isempty(updateDams) %move on if the dams that need updating have been removed also
            continue
        end
        downDam = rmDams(i);
        while gone %search downstream for the next dam that still exists (or the watershed pour point if no downstream dam exists)
            newDownDam = DamIndex(DamIndex(1:end-island,1)==downDam,2);
            if island && k<3 %pass island dam index down to next present dam
                if DamIndex(DamIndex(1:end-island,1)==downDam,3)>0
                    ActiveDamIndex(ismember(DamIndex(1:end-island,1),updateDams),3)=DamIndex(DamIndex(1:end-island,1)==downDam,3);
                    ActiveDamIndex(DamIndex(1:end-island,1)==downDam,3)=0;
                end
            end
            if ismember(newDownDam,rmDams)
                prevDownDam=downDam;
                downDam = newDownDam;
            else
                ActiveDamIndex(ismember(DamIndex(1:end-island,1),updateDams),k)=newDownDam; %update the active dam list
                break
            end
        end
    end
end
ActiveDamIndex(x==0,:)=[];
if island
    %col 3: remove removed dams from col 3. Necessary here because they get
    %shifted to the appropriate downstream location.
    ActiveDamIndex(ismember(ActiveDamIndex(:,3),rmDams),3)=0;
end
%For islands: clear all column four values except for one of the island
%dams just below the island node dam (milford for Penob, 7435). The island
%node can be found at the end, column 1 of DamIndex/ActiveDamIndex.
%calltype, a binary, indicates whether this is an active run, or prep for
%an active run. If it's prep (0), you want to keep flagging all of the island
%dams. If it's active (1), you want only one island dam (adjacent downstream 
%east i.e. third column match) to have the flag.
if island && calltype
    if any(ActiveDamIndex(:,3)>0 & any(any(ActiveDamIndex(ActiveDamIndex(:,3)>0,2:3)>1e9)))
        for i=1:length(ActiveDamIndex)-2
            if ActiveDamIndex(i,3)>0 && any(ActiveDamIndex(i,2:3)>1e9)
                ActiveDamIndex(i,2:3)=[min(ActiveDamIndex(i,2:3)) 0];
            end
        end
%         ActiveDamIndex(ActiveDamIndex(:,3)>0 & any(ActiveDamIndex(ActiveDamIndex(:,3)>0,2:3)>1e9),2:3)=[min(ActiveDamIndex(ActiveDamIndex(:,3)>0 & any(ActiveDamIndex(ActiveDamIndex(:,3)>0,2:3)>1e9),2:3)) 0];
    end
    if sum(ActiveDamIndex(:,3))==0 %sum(ActiveDamIndex(:,4)>1)<1 || 
        ActiveDamIndex(:,4)=0;
    else
        try
            ActiveDamIndex(~ismember(ActiveDamIndex(:,1),ActiveDamIndex(:,3)),4)=0;
            %ActiveDamIndex(ActiveDamIndex(:,1)~=ActiveDamIndex(:,3),4)=0;
        catch
            fprintf('heyo')
        end
    end
end
% here's the old version to correct for islands. Not sure what I meant by it.
% if island
%     IslandDams=find(ActiveDamIndex(:,4)>0);
%     for j=1:length(IslandDams)
%         if ActiveDamIndex(ActiveDamIndex(:,1)==ActiveDamIndex(IslandDams(j),2),4)>0
%             ActiveDamIndex(IslandDams(j),4)=0;
%         end
%     end
% end
