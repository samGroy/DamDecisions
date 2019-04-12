WID = 1000000021;%7487; %
IDs=dlmread(sprintf('%i.txt',WID));
%IDs(IDs==4489)=[];
%Remove watershed IDs that don't match your base watershed ID
%IDs(IDs >= 1e9 && IDs ~= WID)=[];
UpMatrix=zeros(length(IDs),length(IDs));

fprintf('Reading data for basin %i\n',WID)
for i=1:length(IDs)
    try
        d=dlmread(sprintf('%i.txt',IDs(i)));
        %d(ismember(d,4489))=[]; %hack to remove pesky Gilman Falls dam node, necessary so the Penob Island dams don't implode
    catch
        d=IDs(i);%true if watershed is nonexistant or so small that even its own point does not fit into it
    end
    if length(d)>length(IDs)
        continue
    end
    d(d==IDs(i))=0;
    UpMatrix(1:length(d),i)=d;
end

UpCheck=zeros(length(IDs),length(IDs));

fprintf('Removing endless loop errors\n')
for i=1:length(IDs)
    %code to remove endless loops due to barrier wshds of identical size
    %b=IDs(ismember(IDs,UpMatrix(UpMatrix(:,i)>0,i))); %This one seemed convoluted
    b=UpMatrix(UpMatrix(:,i)>0,i);%list IDs of dams upstream of dam i
    if ~isempty(b) %if there is at least one upstream dam (that didn't match IDs(i),already removed from list)
        %col=find(UpMatrix(:,b)==IDs(i));
        a=find(ismember(IDs,b)); %a = indices of the IDs listed in b
        [row,col]=find(UpMatrix(:,a)==IDs(i)); %find index location of the barriers in UpMatrix that cause endless loop: downstream barrier lists upstream barrier as downstream when their watershed polygons are identical, occurs when the two barriers are closer than the resolution of the watershed delineation.
        if ~isempty(col)
            UpMatrix(row,a(col))=0; %remove 'downstream' barrier ID from 'upstream' barrier
            %for k=1:length(col)
             %   UpMatrix(UpMatrix(:,IDs(i))==b(col(k)),IDs(i))=0;
            %end
        end
    end
end

fprintf('Finding adjacent unique upstream pairs\n')
for i=1:length(IDs)
    UpTemp=UpMatrix;
    [~,col]=find(UpMatrix==IDs(i)); %all dams that have i upstream
    col(end+1)=i; %including dam i
    UpTemp(:,col)=[]; %remove them from the temp matrix
    UpCheck(:,i)=~ismember(UpMatrix(:,i),UpTemp(:)); %keep upstream dams if they do not show up for any of the other dams
end

[r,c]=find(UpCheck==1);
for i=1:length(r)
    UpCheck(r(i),c(i))=UpMatrix(r(i),c(i));
end

fprintf('Listing downstream pairs\n')
for i=1:length(IDs)
    ids=UpCheck(UpCheck(:,i)>0,i);
    if isempty(ids)
        continue
    else
        row=find(ismember(IDs(:,1),ids));
        IDs(row,2)=IDs(i);
    end
end
d=1;
dlmwrite('DefaultDownDamList.txt',IDs,'precision','%9.f');
dlmwrite('DefaultUpstreamMatrix.txt',UpMatrix,'precision','%9.f');
fprintf('Completed, but you will need to add island barriers manually\n')
