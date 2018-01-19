
%%% implement x correlation in a c way.
function [out]= myxcorr(in1,in2)
if length(in1)>length(in2)
    pad = length(in1)-length(in2);
    in2 = [in2 zeros(pad,1)];
elseif length(in1)<length(in2)
    pad = length(in2)-length(in1);
    in1 = [in1 zeros(pad,1)];
end

out_len = length(in1);
out = zeros(1,out_len);

for k = 1:out_len
    sum = 0;
    for i = 1:out_len-k
        sum = sum + in1(i+k-1)*in2(i);    
    end
    out(k)=sum;
end
plot(out);