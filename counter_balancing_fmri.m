clear

run = 1:9;
syn = ["O","S"];
cla = ["clear","noisy1","noisy2"];
spr = ["0.75","1.0","1.25"];
ntrial = 9*2*3*3;

cond = strings(ntrial,6);
in = 0;
for r=1:9
    for s=1:2
        for c=1:3
            for p=1:3
                in = in+1;
                cond(in,1:4) = [run(r), syn(s), cla(c), spr(p)];
            end
        end
    end
end

%%

for s=1:2
    str = CreateLatinSquare(9);
    str = str(randperm(9),:);
    in = 0;
    for c=1:3
        for p=1:3
            in=in+1;
            idx = ( cond(:,2)==syn(s) & cond(:,3)==cla(c) & cond(:,4)==spr(p) );
            cond(idx,5) = str(in,:);
        end
    end
end

gen = ["F","M"];
for t=1:ntrial
    cond(t,6) = strcat("00",cond(t,5),"_f",cond(t,2),gen(ceil(2*rand)),"_",cond(t,4));
end


