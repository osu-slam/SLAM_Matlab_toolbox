function str = num2digit(num,digit)
% convert a number into string with a specific number of digit
%
% input arguments:
% num -- input number (double)
% digit  -- n digits (>1)
% 
% output arguments:
% str -- output number (string)

digit = round(digit);
if digit<2
    error('the number of digit is less than 2')
end

if isscalar(num)
    str=num2str(num);
    for d=2:digit
        if num<10^(d-1)
            str = ['0' char(str)];
        end
    end
else
    str = strings(size(num));
    for i=1:size(num,1)
        for j=1:size(num,2)
            str(i,j) = num2str(num(i,j));
            for d=2:digit
                if num(i,j)<10^(d-1)
                    str(i,j) = ['0' char(str(i,j))];
                end
            end
        end
    end
end



end