

syn = ["OO","OS","SO","SS"];
cla = ["clear","noisy"];
tg = ["M", "F"];
fg = ["M", "F"];

cond = strings(96,4);
in = 0;
for rep=1:6
    for s=1:4
        for c=1:2
            for g=1:2
                in = in+1;
                cond(in,1:3) = [syn(s), cla(c), tg(g)];
            end
        end
    end
end

b1 = ["M"; "M"; "F"; "F"];
b2 = [ "F"; "F"; "M"; "M"];
gseq = repmat([b1;b2;b2;b1],3,1);
cond(1:48,4) = gseq;

gseq = repmat([b2;b1;b1;b2],3,1);
cond(49:96,4) = gseq;
