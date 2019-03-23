function out =util_antenna_based_block_sp_ar(u,ant2,param)
param.pos=[0 param.pos];
out.partition = [] ;
length(param.pos)-1
for j =1:length(param.pos)-1
    
blockSZ = max(1, floor(param.pos(j+1)/param.size));
blockMod = mod(param.pos(j+1),param.size);
if (blockMod > (param.size *0.5))
   blockSZ = blockSZ+1;
end

last = sum(param.pos(1:j+1));
if j>1
first = sum(param.pos(1:j))+1;
else
    first =1;
end
ant2__ = ant2(first:last);

u__ = u(first:last);
M = length(ant2__);
%%
t =1;
cmt = 1;
snap_sz =[];
snaps = [];
start = 10;
for  i = start: length(ant2__)
    if ant2__(i) >ant2__(i-1) %new snapshot
        snap_sz(t)= i-cmt;
        cmt =  i;
        snaps(i-1) = t;
        t=t+1;       
    else
        snaps(i-1) =t;
    end
end
snaps(1:start-1) = snaps(start);
snap_sz(t) = M - cmt;

if length(snap_sz)~=snaps(end)
    disp '!1 ERROR SNAPSHOT DELIMITING'
    break;
end
%%
snap_maxu = zeros(size(snap_sz));
snap_maxu(1)= max(abs(u__(start:snap_sz(1))));
for i = 2:length(snap_sz)
    
    indexm = sum(snap_sz(1:(i-1))) + 1;
    indexM = sum(snap_sz(1:i));
    snap_maxu(i)= max(abs(u__(indexm:indexM)));
end

diffu = [];
for i = 2:length(snap_maxu)
    diffu(i)= abs(snap_maxu(i)-snap_maxu(i-1));   
end
diffu(1) =0;

snap_transit_index = find(diffu> 0.1* min(mean(diffu),median(diffu)));

ustair =[];
ustair_sz = [];
t = 1;
for i =1:length(snap_transit_index)
    tend = sum(snap_sz(1:snap_transit_index(i)-1));
    ustair(t:tend) = i; %ustair(t) = -1; ustair(tend) = -1;
    ustair_sz(i) = tend-t+1;
    t = tend +1;
end
ustair_sz(i+1) = M- sum(ustair_sz(1:i));

%%
partition =[];

if param.snapshot ==1
   snapNbr = floor(length(snap_sz)/blockSZ);
   t = 1;
   for i =1:blockSZ
       partition(i) = sum(snap_sz(t:min(t+snapNbr-1,length(snap_sz))));
       fact =floor( partition(i)/param.size);
       fmod =  mod(partition(i),param.size);
       if (fmod > param.size *0.8 ) && (fact ==0)          
           fact =1;
       end
       if (fact ==1)
           if   (fmod < param.size *0.1)
              t = t+snapNbr;
           else
               fact =2;
           end
       end   
       if fact >1
          kk=1;
          limit= 0 ;
          while ((fact>1) && (kk< snapNbr)) || (fact ==1 && fmod > param.size *0.1)
            limit = t+snapNbr-kk;
            partition(i) = sum(snap_sz(t:limit-1));
            fact =floor( partition(i)/param.size);
            fmod = mod(partition(i),param.size);
            kk = kk+1;           
          end
          t = t+snapNbr-kk -1;       
       elseif fact ==0
           kk = 1;       
           while fact ==0
              limit = t+snapNbr+kk;
              partition(i) = sum(snap_sz(t:limit-1));
              fact =floor( partition(i)/param.size);
              fmod =  mod(partition(i),param.size);
               if (fmod > param.size *0.99 ) && (fact ==0)
                  fact = 1;
               else
                  kk = kk+1;
               end
           end
           t = limit;       
       end
   end
   if sum(partition)<param.pos(j+1)
      rest_  = param.pos(j+1) - sum(partition);
      if rest_ < ( mean(partition) *0.8)
        partition(end) = partition(end) +rest_;
      else 
        partition = [partition rest_];
      end
   end  
   
else
    max_nblocks = length(ustair_sz);
%     if blockSZ >= max_nblocks
%        partition = ustair_sz;
%     else
     blockSZ = min(max_nblocks,blockSZ);
        
        t =1;
        i =0;
        tend =0;
        while i < blockSZ+1 && t < max_nblocks
            i=i+1;
            etemp  = sum(ustair_sz(t:t));
            stepv = 0;
            
            
            while etemp < 0.99*param.size  &&  t+stepv<max_nblocks          
                stepv = stepv +1;                
                etemp  = sum(ustair_sz(t: t+stepv));
                tend = t+stepv;
            end   
            
            partition = [partition etemp];
            
            t = t+stepv+1;
            rest_ = M - sum(partition);
            rest_F = floor(rest_/param.size);
            rest_M = mod(rest_,param.size);
            
            if (rest_F ==0 && rest_M < param.size *0.5)
                
                i =  blockSZ+1;
            end
        
        end
        
        
        sz_tmp =length(partition);
        rest_ = M - sum(partition);
        if rest_>0
           if rest_M > param.size * 0.75
               partition = [partition rest_] ;
           else
               partition(sz_tmp) = partition(sz_tmp)+ rest_ ; 
           end
          
        end
        
%     end
    
    
    
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
