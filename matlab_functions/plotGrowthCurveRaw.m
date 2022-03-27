function plotGrowthCurveRaw(filename, varargin)

[~,basename,ext] = fileparts(filename);

%%
p = inputParser;
addRequired(p,'filename',@ischar);
addParameter(p,'samples','',@ischar);
addParameter(p,'Save','');
addParameter(p,'Exclude','');

parse(p,filename,varargin{:});
samples = p.Results.samples;
saveFormat = p.Results.Save;
excluded = p.Results.Exclude;
%%
switch lower(ext)
    case '.asc'
        % read data from ASC file (CSV)
        data = readtable(p.Results.filename,'Filetype','text','ReadVariableNames',true,...
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
    case '.txt'
        % read data from TXT file (TAB)
        data = readtable(p.Results.filename,'ReadVariableNames',true);
        data = table2struct(data,'ToScalar',true);
        data.sample = string(data.sample);
        data.position = string(data.replicate);
    otherwise
        error('File format not recognized');
end
%%

ind = contains(data.sample,samples);
data = data(ind,:);

wells = unique(data.position);
wells = wells(~ismember(wells,excluded));
newdata = [];

for k = 1:numel(wells)
    ind = strcmp(data.position,wells{k});
    subdata = data(ind,:);
    [~,order] = sort(subdata.time);
    subdata = subdata(order(2:end-1),:);
    
    % rate
    dOD = diff(subdata.OD)./diff(subdata.time);
    subdata.dOD = ([dOD(1); dOD] + [dOD; dOD(end)])/2 ./subdata.OD;
    ii = find(subdata.OD<=0,1,'last');
    subdata.dOD(1:ii) = NaN;
    
    newdata = [newdata; subdata];
end

%%
samples = unique(newdata.sample);
nb = regexp(samples,'[0-9]+','match');
if numel(nb)>1
    [~,nb] = sort(str2double(cellfun(@(x) strjoin(x,''), nb, 'UniformOutput', false)));
    samples = samples(nb);
end
%%

for k = 1:numel(samples)
    
    ind = strcmp(newdata.sample,samples{k});
    data = newdata(ind,:);
    
    [pos,~,pc] = unique(data.position);
    c = lines(numel(pos));
    
    h1 = figure('name',strcat(basename,'_',samples{k}),'position',[200 200 800 800]);
    
    
    ax1 = subplot(2,2,1); hold on;
    
    r =zeros(numel(pos),1);
    for j = 1:numel(pos)
        r(j) = scatter(ax1,data.time(pc==j),data.OD(pc==j),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
    end
    
    axis square; axis tight;
    xlabel('Time (min)');
    ylabel('OD_{600}');
    title(ax1,strcat(basename,'_',samples{k}),'Interpreter','none');
    legend(ax1,r,cellstr(pos),'Location','NorthWest','Fontsize',10);
    ax1.FontSize = 12;
    
    ax2 = subplot(2,2,2); hold on;
    
    for j = 1:numel(pos)
        scatter(ax2,data.time(pc==j & data.OD>0),data.OD(pc==j & data.OD>0),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
    end
    ax2.XLim = ax1.XLim;
    axis square; axis tight;
    xlabel('Time (min)');
    ylabel('OD_{600}');
    ax2.FontSize = 12;
    ax2.YScale = 'log';
    
    
    ax3 = subplot(2,2,3); hold on;
    
    for j = 1:numel(pos)
        scatter(ax3,data.time(pc==j & data.OD>0),data.dOD(pc==j & data.OD>0),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
    end
    
    ax3.XLim = ax1.XLim;
    axis square; axis tight;
    xlabel('Time (min)');
    ylabel('Growth rate (min^{-1})');
    ax3.FontSize = 12;
    plot(ax3,ax3.XLim, [0 0],'k');
    
    
    ax4 = subplot(2,2,4); hold on;
    
    for j = 1:numel(pos)
        scatter(ax4,data.OD(pc==j & data.OD>0),data.dOD(pc==j & data.OD>0),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
        
    end
    axis square; axis tight;
    xlabel('OD_{600}');
    ylabel('Growth rate (min^{-1})');
    ax4.FontSize = 12;
    ax4.XScale = 'log';
    plot(ax4,ax4.XLim, [0 0],'k');
    
    if ismember('svg',saveFormat)
    set(h1,'PaperPositionMode','auto');
    print(h1,sprintf('%s_%s.svg',basename,samples{k}),'-dsvg','-r300','-noui');
    end


    if ismember('png',saveFormat)
        set(h1,'PaperPositionMode','auto');
        print(h1,sprintf('%s_%s.png',basename,samples{k}),'-dpng','-r300','-noui');
        crop(sprintf('%s_%s.png',basename,samples{k}));
    end
    
    
end

%%

end