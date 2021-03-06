function nets = DPPFwkshp_fitfn_split(x,DamIndex,OVs,VarNames,OVcount)
%This code is like 'DPPFwkshp_fitfn' but splits the fish habitat results
%per dam. Pretty sure it's just viciously hacked as I ran out of time.
x=floor(x);
if OVs{1,end}(8)>=4
    t=1;
elseif OVs{1,end}(8)==3 && size(OVs{strcmp('power',VarNames)},2)==1
    x(x==2)=3;
    t=1;
elseif OVs{1,end}(8)==3 && size(OVs{strcmp('power',VarNames)},2)==2
    t=0;
else
    t=0;
end

%potential hydropower development sites: under status quo wouldn't exist,
%so they can only be triggered if x is 2 or 4
%make sure their hydropower is in 2nd column, installation cost is in 2nd
%column
if size(OVs{strcmp('power',VarNames)},2)==2
    x(DamIndex(1:length(x),1)'>=9e8 & ismember(x,[1 3]))=0;
end
% %cost hack
% x(x>0)=3;
out = zeros(size(x,1),OVcount);

%NSGA-II has generated a new population of x values, now revise the
%DamIndex to these changes


%revise DamData (aka ObjectiveVars) ased on new x
% ActiveObjectiveVars=ObjectiveVars;
% for i=1:size(ActiveObjectiveVars,2)
%     ActiveObjectiveVars{1,i}(x==0,:)=[];
% end

%One thing is checked and two things are revised here:
%1. Check if the watershed has island dams
%2. Get the new active dam index, the lenth of which differs if the island
%dams exist, and in all cases the index is longer than x by at least one
%because the watershed ID (at the bottom) is excluded from selection in
%NSGA-II as it is a constant feature.
%3. Simple sum functions, set the proper datarange before summing by
%removing the watershed and island dam data at the end
% if size(DamIndex,2)>2
datarange=length(x);
%x(x(DamIndex(1:length(x),1)>9e8)==1)=0; %constructing new dams: if x==1, make x==0 because for status quo, the dam doesn't exist
x = [x ones(1,length(DamIndex)-length(x))];
ActiveDamIndex=DPPF_idx(x,DamIndex,1);
% else
%     x = [x 1];
%     ActiveDamIndex=DPPF_idx(x,DamIndex);
%     datarange=size(DamIndex,1)-1;
% end

ct=0; %count off OVs
if ismember('power',VarNames)
    if size(OVs{strcmp('power',VarNames)},2)==1
        out(7) = -1*sum(OVs{strcmp('power',VarNames)}(x(1:datarange)>0)); %hydropower
    else
        out(7) = -1*(sum(OVs{strcmp('power',VarNames)}(x(1:datarange)>0,1))+...
            sum(OVs{strcmp('power',VarNames)}(ismember(x(1:datarange),[2,4]),2))); %hydropower refurb and new
    end
end

if ismember('storage',VarNames)
    out(3) = -1*sum(OVs{strcmp('storage',VarNames)}(x(1:datarange)>0)); %reservoir volume
end

if ismember('property',VarNames)
    out(6) = -1*sum(OVs{strcmp('property',VarNames)}(x(1:datarange)==1)); %disturbed properties
end

if ismember('cost',VarNames)
    %keep costs
    out(4)=out(4)+sum(OVs{strcmp('cost',VarNames)}(x(1:datarange)==1,1));
    %remove costs
    out(4)=out(4)+sum(OVs{strcmp('cost',VarNames)}(x(1:datarange)==0,5));
    %improve hydro costs
    out(4)=out(4)+sum(OVs{strcmp('cost',VarNames)}(x(1:datarange)==2,2));
    %improve fish costs
    out(4)=out(4)+sum(OVs{strcmp('cost',VarNames)}(x(1:datarange)==3,3));
    %improve both costs
    out(4)=out(4)+sum(OVs{strcmp('cost',VarNames)}(x(1:datarange)==4,4));
end

if ismember('hazard',VarNames)
    %function 2: return the number of high risk dams
%     out(OVindex(strcmp('hazard',VarNames)>0,1):OVindex(strcmp('hazard',VarNames)>0,2)) = sum(ObjectiveVars{strcmp('hazard',VarNames)}(x(1:datarange)==1)==3); %dam hazard level    
    out(5) = -1*sum(OVs{strcmp('hazard',VarNames)}(x(1:datarange)==0)); %high hazard dams removed  
end

if ismember('CO2',VarNames)
    out(8) = -1*sum(OVs{strcmp('CO2',VarNames)}(x(1:datarange)>0)); %dams that remain offset CO2
end

if ismember('indigenous',VarNames)
    %keep
    out(9)=out(9)+sum(OVs{strcmp('indigenous',VarNames)}(x(1:datarange)==1,1));
    %remove
    out(9)=out(9)+sum(OVs{strcmp('indigenous',VarNames)}(x(1:datarange)==0,5));
    %improve hydro
    out(9)=out(9)+sum(OVs{strcmp('indigenous',VarNames)}(x(1:datarange)==2,2));
    %improve fish
    out(9)=out(9)+sum(OVs{strcmp('indigenous',VarNames)}(x(1:datarange)==3,3));
    %improve both
    out(9)=out(9)+sum(OVs{strcmp('indigenous',VarNames)}(x(1:datarange)==4,4));
    out(9)=-out(9);
end

if ismember('historic',VarNames)
    %keep
    out(10)=out(10)+sum(OVs{strcmp('historic',VarNames)}(x(1:datarange)==1,1));
    %remove
    out(10)=out(10)+sum(OVs{strcmp('historic',VarNames)}(x(1:datarange)==0,5));
    %improve hydro
    out(10)=out(10)+sum(OVs{strcmp('historic',VarNames)}(x(1:datarange)==2,2));
    %improve fish
    out(10)=out(10)+sum(OVs{strcmp('historic',VarNames)}(x(1:datarange)==3,3));
    %improve both
    out(10)=out(10)+sum(OVs{strcmp('historic',VarNames)}(x(1:datarange)==4,4));
    out(10)=-out(10);
end

if ismember('town',VarNames)
    %keep
    out(11)=out(11)+sum(OVs{strcmp('town',VarNames)}(x(1:datarange)==1,1));
    %remove
    out(11)=out(11)+sum(OVs{strcmp('town',VarNames)}(x(1:datarange)==0,5));
    %improve hydro
    out(11)=out(11)+sum(OVs{strcmp('town',VarNames)}(x(1:datarange)==2,2));
    %improve fish
    out(11)=out(11)+sum(OVs{strcmp('town',VarNames)}(x(1:datarange)==3,3));
    %improve both
    out(11)=out(11)+sum(OVs{strcmp('town',VarNames)}(x(1:datarange)==4,4));
    out(11)=-out(11);
end

if ismember('aesthetic',VarNames)
    %keep
    out(12)=out(12)+sum(OVs{strcmp('aesthetic',VarNames)}(x(1:datarange)==1,1));
    %remove
    out(12)=out(12)+sum(OVs{strcmp('aesthetic',VarNames)}(x(1:datarange)==0,5));
    %improve hydro
    out(12)=out(12)+sum(OVs{strcmp('aesthetic',VarNames)}(x(1:datarange)==2,2));
    %improve fish
    out(12)=out(12)+sum(OVs{strcmp('aesthetic',VarNames)}(x(1:datarange)==3,3));
    %improve both
    out(12)=out(12)+sum(OVs{strcmp('aesthetic',VarNames)}(x(1:datarange)==4,4));
    out(12)=-out(12);
end

%Watershed-dependent variables require their own functions to (for example) subtract
%between dam watersheds and generate the product of fish passage
%probability
%trigger if an objective variable that requires network calculation is included 
nvars=[];
dvars=[];
if ismember('passage',VarNames)
    if ~t %all(x<2) && 
        dvars=[dvars OVs{strcmp('passage',VarNames)>0}(x>0,:)];
    else
        if size(OVs{strcmp('passage',VarNames)},2)==2
            OVs{strcmp('passage',VarNames)>0}(x>2,1)=OVs{strcmp('passage',VarNames)>0}(x>2,2);
            dvars=[dvars OVs{strcmp('passage',VarNames)>0}(x>0,1)];
        else
            OVs{strcmp('passage',VarNames)>0}(x>2,1:2)=OVs{strcmp('passage',VarNames)>0}(x>2,3:4);
            dvars=[dvars OVs{strcmp('passage',VarNames)>0}(x>0,1:2)];
        end
    end
end
if ismember('alewifeBlueback',VarNames)
    nvars=[nvars OVs{strcmp('alewifeBlueback',VarNames)>0}(x>0,:)];
    dvars=[dvars nvars(:,end)];
    nvars(:,end)=[];
end
if ismember('shad',VarNames)
    nvars=[nvars OVs{strcmp('shad',VarNames)>0}(x>0,:)];
end
if ismember('salmon',VarNames)
    nvars=[nvars OVs{strcmp('salmon',VarNames)>0}(x>0,:)];
end

r=0;
if ismember('riverCanoe',VarNames)
    r=r+1;
    nvars=[nvars OVs{strcmp('riverCanoe',VarNames)>0}(x>0,:)];
    dvars=[dvars nvars(:,end-1:end)];
    nvars(:,end-1:end)=[];
end
if ismember('riverRaft',VarNames)
    r=r+1;
    nvars=[nvars OVs{strcmp('riverRaft',VarNames)>0}(x>0,:)];
    dvars=[dvars nvars(:,end-1:end)];
    nvars(:,end-1:end)=[];
end

if sum(any(nvars))
    nets = ...
        -1*DPPF_netOV_split(ActiveDamIndex,nvars,dvars,VarNames,DamIndex);
        %-1*DPPF_netOV(ActiveDamIndex,nvars,dvars,VarNames);
    j=1;
    if any(OVs{strcmp('ag',VarNames)>0}(1:3))
        for i=1:3
            if OVs{strcmp('ag',VarNames)>0}(i)
                if OVs{strcmp('ag',VarNames)>0}(7)
                    out(1)=out(1)+nets(j).*OVs{strcmp('ag',VarNames)>0}(i);
                else
                    ct=ct+1;
                    out(ct)=nets(j).*OVs{strcmp('ag',VarNames)>0}(i);
                end
                j=j+1;
            end
        end
    end
end

if r && OVs{strcmp('ag',VarNames)}(7)>1
    %sum of total non-dam watershed scale flatwater area and the
    %added area provided by present dams (x==1). Total non-dam lake area
    %housed in head node for encompassing watershed (datarange+1).
    ct=max(ct)+1;
    out(ct) = nets(end);
    if r==2
        out(ct) = out(ct) + nets(end-1);
    end
    if ismember('lakeRecreation',VarNames) 
        if OVs{strcmp('ag',VarNames)}(7)>2
            out(ct) = out(ct) -1*(sum(OVs{strcmp('lakeRecreation',VarNames)}(datarange+1:end,1))...
                -sum(OVs{strcmp('lakeRecreation',VarNames)}(x(1:datarange)==0,1)));
        else
            ct=ct+1;
            out(ct) = -1*(sum(OVs{strcmp('lakeRecreation',VarNames)}(datarange+1:end,1))...
            -sum(OVs{strcmp('lakeRecreation',VarNames)}(x(1:datarange)==0,1)));
        end
    end
else
    if r
        out(2)=nets(end);%/OVs{strcmp('ag',VarNames)}(6)/OVs{strcmp('ag',VarNames)}(5);
        if r==2
            ct=ct+1;
            out(ct)=nets(end-1);%/OVs{strcmp('ag',VarNames)}(6)/OVs{strcmp('ag',VarNames)}(5);
        end
    end
    if ismember('lakeRecreation',VarNames)
        ct=ct+1;
        out(ct) = -1*(sum(OVs{strcmp('lakeRecreation',VarNames)}(datarange+1:end,1))...
            -sum(OVs{strcmp('lakeRecreation',VarNames)}(x(1:datarange)==0,1)));%...
            %/OVs{strcmp('ag',VarNames)}(6)/OVs{strcmp('ag',VarNames)}(4);
    end
end