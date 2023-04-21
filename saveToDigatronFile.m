function saveToDigatronFile(timeDuration,signal,nvp)

arguments
    timeDuration {mustBeVector}
    signal {mustBeVector}
    nvp.signal = "I"
    nvp.saveFileName = "fileName.txt"
end

fileID = fopen(nvp.saveFileName,'w');
for tt = 1:numel(timeDuration)
    if strcmp(nvp.signal,"I")
        fprintf(fileID,"%1.0f sec;%0.4f;;;\r\n",timeDuration(tt),signal(tt));
    elseif strcmp(nvp.signal,"P")
        fprintf(fileID,"%1.0f sec;;%0.4f;;\r\n",timeDuration(tt),signal(tt));
    end
end
fclose(fileID);

end