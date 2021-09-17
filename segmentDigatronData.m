function dataSeg = segmentDigatronData(data,idx)
numFields = fields(data);

for nF = 1:length(numFields)
    tempSig = data.(numFields{nF});
    tempSig = tempSig(idx(1):idx(2));
    if contains(numFields{nF},{'ProgTime'})
        tempSig = tempSig - tempSig(1);
    end
    dataSeg.(numFields{nF}) = tempSig;
end
end