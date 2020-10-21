function out =util_time_based_block_sp_ar(u,time_,param)
param.pos=[0 param.pos];
out.partition = [] ;
length(param.pos)-1
for j =1:length(param.pos)-1
    
blockSZ = max(1, floor(param.pos(j+1)/param.size));
blockMod = mod(param.pos(j+1),param.size);
if (blockMod > (param.size *0.5)) && (param.size<param.pos(j+1))
   blockSZ = blockSZ+1;
end

last = sum(param.pos(1:j+1));
if j>1
first = sum(param.pos(1:j))+1;
else
    first =1;
end
time = time_(first:last);

u__ = u(first:last);
M = length(time);
%%
t =1;
cmt = 1;
snap_sz =[];
snaps = [];
start = 2;
for  i = start: length(time)
    if time(i) ~=time(i-1) %new snapshot
        snap_sz(t)= i-cmt;
        cmt =  i;
        snaps(i-1) = t;
        t=t+1;       
    else
        snaps(i-1) =t;
    end
end
snaps=[snaps snaps(end)];
snaps(1:start-1) = snaps(start);
snap_sz(t) = M - cmt+1;

if length(snap_sz)~=snaps(end)
    disp '!1 ERROR SNAPSHOT DELIMITING'
    break;
end
%%

for i=1:length(snap_sz)-1
    
    snapEndTime=time(sum(snap_sz(1:i)));
    snapStartTime=time(sum(snap_sz(1:i))+1);
    diffSnapTime(i)  = snapStartTime -snapEndTime;
end
[~,startTrack]=find(abs(diffSnapTime)>3*median(abs(diffSnapTime)));

for i=1:length(startTrack)
    trackStartPos(i)=sum(snap_sz(1:startTrack(i)))+1;
end


%%
partition =[];

    max_nblocks = length(startTrack)+2;
    blockSZ = min(max_nblocks,blockSZ);
    t =1;
    i =0;
    snapStart=1;
    snapLast=startTrack(t);
    if (blockSZ < max_nblocks)
       while ( i < blockSZ+1 ) &&(t<max_nblocks+1) &&( snapLast-1<length(snap_sz))
           
          tend=sum(snap_sz(snapStart:snapLast));
         
          if tend>param.size
             i=i+1;
             partition=[partition tend];
             
             snapStart=snapLast+1;
             if t> length(startTrack)-1
                snapLast= length(snap_sz);
             else
                 snapLast=startTrack(t+1);
                 t=t+1;
             end
             
          elseif snapLast ==length(snap_sz)
                    if tend>param.size *0.5
                       partition=[partition tend];
                    else
                        partition(end)=partition(end)+tend;
                    end
                    i=inf;                 
           
          elseif  (t<length(startTrack))
                   snapLast= startTrack(t+1);
                   t=t+1;
                   
          elseif (t==length(startTrack))
              snapLast=snapLast+1;
          end
             
       end
       
    else
        
       partition(1)=sum(snap_sz(1:startTrack(1)));
       for i=2:length(startTrack)         
           partition(i)=sum(snap_sz(startTrack(i-1)+1:startTrack(i)))
       end
       partition(length(startTrack)+1)=sum(snap_sz(startTrack(end)+1:end));
       
    end
             
             
        


if sum(partition)~=param.pos(j+1)
    disp '!3 ERROR IN BLOCK SPLITTING'
    break;
end

out.blockNumber(j) = length(partition);
out.partition =[out.partition partition];

end
disp(out.partition)
end
