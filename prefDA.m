function pref=prefDA(dat,x)
%based on AHP input from users, generate global weighted preferences based
%on raw preferences and the frequency in which each decision alternative is
%used in the MOGA output
pref=zeros(size(x,1),size(dat,2)); %preferences for each criteria, for each scenario
freq=zeros(1,size(dat,1));%-1); %frequency of each decision alternative, top row is global so subtract one from length
da=[0 3 2 4 1];
%loop through scenarios to find frequency of each decision alternative
for i=1:size(pref,1)
    %loop through each decision alternative to find their frequency
    for j=1:length(freq)
        freq(da(j)+1)=sum(x(i,:)==da(j))/size(x,2); %fraction decisions for this specific alternative, all decisions within the scenario.
    end
    pref(i,:)=sum(freq'.*dat,1);%.*dat(1,:); %sum across decision alternatives. multiply frequency by local weights. 2019 test didn't incorporate global weights, so do not include here. Multiple result by global weight, one per criterion.
    pref(i,:)=pref(i,:)./sum(pref(i,:));%level out to unity.
end
