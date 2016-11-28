function output=error_ratio(input,label)
output=sum(sum(input~=label))/size(label,2);
end