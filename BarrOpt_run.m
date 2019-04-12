function BarrOpt_run(wshd,ovs,pop,x)
%BarrOpt_run.m

if ispc
    basedir = 'D:\FoD\PPF';
else
    basedir = '//home/sroy/FoD/PPF';
end

for j=1:size(ovs,1)
    ovars=ovs(j,:);
    [BarrIndex, v] = BarrOpt_prep(wshd,basedir,'crh1','zero',ovars);
    if nargin<4
        if size(BarrIndex,2)>2
            island=1;
        else
            island=0;
        end
        x=ones(2,length(BarrIndex(BarrIndex(:,1)<1e9,1))-island);
        x(2,:)=0;
    elseif iscell(x)
        load(x{1},x{2});
    end
    
    fprintf('run 1 of 3\n')
    [x, f, ~, ~, ~, ~, ~, ~]=BarrOpt(BarrIndex,pop,x,v);
    fprintf('run 2 of 3\n')
    [x, f, ~, ~, ~, ~, ~, ~]=BarrOpt(BarrIndex,pop,x,v);
    fprintf('run 3 of 3\n')
    [x, f, xu, fu, exitflag, output, population, score]=BarrOpt(BarrIndex,pop,x,v);
    
    tm=datestr(now);
    tm=strrep(tm,':','-');
    tm=strrep(tm,' ','_');
    o=[ovars{1}];
    for i=2:length(ovars)
        o=[o '_' ovars{i}];
    end
    if iscell(wshd)
        w=[wshd{1}];
        for i=2:length(wshd)
            w=[w '_' wshd{i}];
        end
    else
        w=wshd;
    end
    save(sprintf('%s/BarrOut/%s_%s_%s_%i.mat',basedir,w,o,tm,pop));
end