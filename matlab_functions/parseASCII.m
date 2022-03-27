function data = parseASCII(filename, varargin)
        data = readtable(filename,'Filetype','text','ReadVariableNames',true,...
            'ReadRowNames',true,'MultipleDelimsAsOne',false);
        ind = find(strcmp(data.Properties.VariableNames,'RawData'));
        samplingtime = data(1,ind:end-1);
        samplingtime = str2double(strrep(samplingtime.Variables,'s',''))/60;
        temperature = data(2,ind:end-1);
        temperature = str2double(regexprep(temperature.Variables,'[^\d\.]',''));
        rawdata = data(3:end,ind:end-1);
        rawdata = cellfun(@str2double,rawdata.Variables);

        metadata = data(3:end,1:ind-1);
        sample = string(repmat(metadata.Layout,1,size(rawdata,2)));
        position = string(repmat(metadata.Properties.RowNames,1,size(rawdata,2)));
        samplingtime = repmat(samplingtime,size(rawdata,1),1);
        temperature = repmat(temperature,size(rawdata,1),1);
        data = struct('sample',sample(:),...%'replicate',replicate(:),...
            'time',samplingtime(:),'OD',rawdata(:),'temperature',...
            temperature(:),'position', position(:));
        data = struct2table(data);
        writetable(data, strcat(filename, ".csv"));
end