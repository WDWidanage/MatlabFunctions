
function [arrayValue, arrayPos] = findClosest(values, array)

for ii = 1:length(values)
    d = abs(array(:) - values(ii));
    [~,idxClosest] = min(d);
    arrayValue(ii) = array(idxClosest);
    arrayPos(ii) = idxClosest;
end

end