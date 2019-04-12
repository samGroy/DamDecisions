%[ba, va] = BarrOpt_prep('all','D:\FoD\PPF','crh1','zero',ovs);
%v=va;
%norm function
%va=v;
scalecode='TOWNS';
amt=0.0522e6;%TOWNS 18.18e6;%HUC08

[ba, va] = BarrOpt_prep('all','D:\FoD\PPF','crh1','zero',ovs);
xr=ones(1,length(ba(:,1))-6);

d=dir([scalecode '-' '*' '.mat']);

for k=1:length(d)
    load([d(k).folder '/' d(k).name]);
    v=va;
    xa=ones(size(x,1),length(ba(:,1)));
    [tf,idx]=ismember(BarrIndex(:,1),ba(:,1));
    if size(BarrIndex,2)>2
        sbtr=6;
    else
        sbtr=5;
    end
    xa(:,idx(1:end-sbtr))=x;
    out=zeros(size(x,1),2);
    
    for i=1:length(out)
        out(i,:)=-BarrOpt_fitfn(xa(i,:),ba,v);
    end
    
    df=(out(:,2)+amt).^2;
    a=find(df==min(df));
    xr(1,idx(1:end-sbtr))=x(a(1),:);
end