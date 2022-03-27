function fitGrowthCurveSTAN(filename, varargin)

[~,baseName,ext] = fileparts(filename);

%%
p = inputParser;
defaultChain  = 4;
defaultWarmup = 500;
defaultThin   = 5;
defaultIter   = ceil(500/defaultChain*defaultThin);
defaultModel  = 'generalized_logistic_growth.stan';
defaultOut    = strcat('_',baseName,'_temp.csv');
defaultSample = '';

addRequired(p,'filename',@ischar);
addParameter(p,'chain',defaultChain,@isnumeric);
addParameter(p,'warmup',defaultWarmup,@isnumeric);
addParameter(p,'thin',defaultThin,@isnumeric);
addParameter(p,'iter',defaultIter,@isnumeric);
addParameter(p,'modelfile',defaultModel,@ischar);
addParameter(p,'outfile',defaultOut,@ischar);
addParameter(p,'sample',defaultSample);

parse(p,filename,varargin{:});

%%
delete(strcat(strrep(p.Results.modelfile,'.stan',''),'*'));
%%
pathFun = which('fitGrowthCurveSTAN');
[pathFun, ~, ~] = fileparts(pathFun);
modelFolder = fullfile(pathFun(1:end-16),'models');
modelPath   = fullfile(modelFolder, p.Results.modelfile);
copyfile(modelPath);
modelPath = fullfile(pwd, p.Results.modelfile);
workdir = pwd;
if ispc
    usr = getenv('USERNAME');
    curr_usr = upper(strrep(usr,' ',''));
    if length(curr_usr)>6
        curr_usr = strcat(curr_usr(1:6),'~1');
        modelPath = strrep(modelPath,usr,curr_usr);
        workdir = strrep(pwd,usr,curr_usr);
    end
end

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
        replicate = cellfun(@str2double,repmat(regexprep(metadata.ReplicateInfo,'/.*',''),1,size(rawdata,2)));
        samplingtime = repmat(samplingtime,size(rawdata,1),1);
        temperature = repmat(temperature,size(rawdata,1),1);
        data = struct('sample',sample(:),'replicate',replicate(:),...
            'time',samplingtime(:),'od',rawdata(:),'temperature',...
            temperature(:),'position', position);
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
if isempty(p.Results.sample)
    samples = unique(data.sample);
else
    samples = intersect(cellstr(p.Results.sample),unique(data.sample));
end

samples = samples(~contains(samples,'bl','IgnoreCase',true));

for k = 1:numel(samples)
    
    ind = strcmp(data.sample,samples{k});
    
    dataOD.replicate = data.replicate(ind);
    dataOD.time = data.time(ind);
    dataOD.od = data.od(ind);
    dataOD.temperature = data.temperature(ind);
    dataOD.position = data.position(ind);
    
    dataOD_unf = dataOD;
    
    f = fit(dataOD.time,dataOD.od,'smoothingspline');
    tstart = max(dataOD.time(f(dataOD.time) == min(f(dataOD.time))));
    tend = min(dataOD.time(f(dataOD.time) == max(f(dataOD.time))));
    [r, tmax] = max(diff(log2(f(tstart:tend))));
    tmax = tstart+tmax;
    ind = dataOD.time>=tstart & dataOD.time<=tend & dataOD.time>min(dataOD.time) & dataOD.time<max(dataOD.time);
    
    if sum(ind)<5 || f(tend)/f(tstart) < 1.5
        continue;
    end
    
    dataOD.replicate = dataOD.replicate(ind);
    dataOD.time = dataOD.time(ind);
    dataOD.od = dataOD.od(ind);
    dataOD.temperature = dataOD.temperature(ind);
    dataOD.position = dataOD.position(ind);
    
    [~,~,ic] = unique(dataOD.replicate);
    dataOD.replicate = uint8(ic);
    dataOD.N = uint32(numel(dataOD.od));
    dataOD.M = uint32(numel(unique(dataOD.replicate)));
    samplelabel = unique(dataOD.position,'stable');
    
    % Quick estimate parameters
    ft = fittype('A+(1+E*t)*(K-A)/(1+2^((Q-t)/B))','coefficients',{'A','E','K','B','Q'},'independent','t');
    opts = fitoptions(ft);
    opts.StartPoint = [min(dataOD.od) 0 max(dataOD.od) 1/r/2 tmax];
    opts.Lower = [0 0 0 0 0];
    opts.Upper = [Inf Inf Inf Inf Inf];
    out = fit(dataOD.time,dataOD.od,ft,opts);
    %     figure,
    %     plot(dataOD.time,dataOD.od,'o'); hold on; plot(out)
    %     plot(dataOD.time,dataOD.od./(1+out.E.*dataOD.time),'o');
    dataOD.eA = out.A; dataOD.eK = out.K; dataOD.eB = out.B; dataOD.eQ = out.Q;
    if strcmpi(ext,'.asc')
        dataOD.od = (dataOD.od-out.A)./(1+out.E.*dataOD.time)+out.A;
        dataOD_unf.od = (dataOD_unf.od-out.A)./(1+out.E.*dataOD_unf.time)+out.A;
    end
    %     plot(dataOD.time,dataOD.od,'o');
    %     drawnow;
    %          set(gca,'YScale','log');
    %%
    fprintf('Fitting model for %s...\n',samples{k});
    fitOD = stan('file',modelPath,'data',rmfield(dataOD,{'position','temperature'}),'verbose',false,...
        'iter',p.Results.iter,'chains',p.Results.chain,...
        'warmup',p.Results.warmup,'thin',p.Results.thin,...
        'sample_file',strcat(samples{k},p.Results.outfile),...
        'working_dir',workdir);
    
    fitOD.block();
    out = fitOD.extract('permuted',true);
    fprintf('done.\n');
    %% plot data
    c = lines(max(3,dataOD.M));
    t = mean(out.time_pred);
    l = cell(dataOD.M+1,1);
    
    screensize = get( groot, 'Screensize' );
    h = figure('Name',samples{k},'Position',[round((screensize(3)-1500)/2)...
        round((screensize(4)-500)/2) 1500 500]);
    
    subplot(2,5,[1 2 6 7]), hold on,
    ax = gca;
    if sum(~isnan(dataOD.temperature))>0
        yyaxis right;
        g0 = plot(dataOD.time, dataOD.temperature,'*','MarkerEdgeColor',[0.6 0.6 0.6]);
        ylabel('^{\circ}C');
        ylim([floor(min(dataOD.temperature)-2.5) ceil(max(dataOD.temperature)+2.5)]);
        ax.YColor = [0.6 0.6 0.6];
%         l{1} = 'Temperature';
        yyaxis left;
    end
    
    for i = 1:dataOD.M
        ind = dataOD.replicate == i;
        g(2*i-1) = plot(dataOD_unf.time(ind), dataOD_unf.od(ind),'o','MarkerFaceColor',c(i,:),...
            'MarkeredgeColor','none','MarkerSize',8);
        g(2*i) = plot(t,mean(out.OD_pred(:,:,i)),'-','linewidth',2,'color',c(i,:));
        l{i} = sprintf('Sample %s',samplelabel(i));
    end
    g(2*i+1) = plot(t,mean(out.OD_pred_all),'linewidth',3,'color','k');
    l{end} = 'Fit';
    legend(g([1:2:2*dataOD.M 2*dataOD.M+1]),l,'location','NorthWest','FontSize',10);
    
    xlabel('Time (min)','FontSize',12);
    ylabel('OD_{600}','FontSize',12);
    
    ax.FontSize = 10;
    ax.YScale = 'log';
    ax.YLim = [min(dataOD.od) max(dataOD.od)];
    ax.YColor = [0 0 0];
    title(strcat(baseName,': ',samples{k}),'Interpreter','none');
    
    ax.PlotBoxAspectRatio(1) = ax.PlotBoxAspectRatio(2);
    %%
    % plot parameters
    n = double(dataOD.M+1);
    xl = ['All'; samplelabel];
    %
    subplot(2,5,3), hold on,
    A = [out.A out.A_rep];
    %     mode_A = mode(round(A,3,'significant'));
    mean_A = mean(A,1);
    plot(1:size(A,2), mean_A,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    plot(1:size(A,2), quantile(A,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(A,2), quantile(A,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(A,2), median(A,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    %     plot(1:size(A,2), mode_A,'*','MarkerEdgeColor',c(1,:));
    ax = gca;
    ax.XTick = 1:n;
    ax.XTickLabel = xl;
    ax.XLim = [0.5 n+0.5];
    ax.YLabel.String = 'OD_{600}';
    ax.Title.String = sprintf('A (mean=%G OD)',round(mean_A(1),3,'significant'));
    ax.FontSize = 8;
    %
    subplot(2,5,4), hold on,
    K = [out.K out.K_rep];
    %     mode_K = mode(round(K,3,'significant'));
    mean_K = mean(K,1);
    plot(1:size(K,2), mean_K,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    plot(1:size(K,2), quantile(K,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(K,2), quantile(K,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(K,2), median(K,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    %     plot(1:size(K,2), mode_K,'*','MarkerEdgeColor',c(1,:));
    ax = gca;
    ax.XTick = 1:n;
    ax.XTickLabel = xl;
    ax.XLim = [0.5 n+0.5];
    ax.YLabel.String = 'OD_{600}';
    ax.Title.String = sprintf('K (mean=%G OD)',round(mean_K(1),3,'significant'));
    ax.FontSize = 8;
    %
    subplot(2,5,5), hold on,
    B = exp([out.B out.B_rep]);
    %     mode_B = mode(round(B,3,'significant'));
    mean_B = mean(B,1);
    plot(1:size(B,2), mean_B,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    plot(1:size(B,2), quantile(B,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(B,2), quantile(B,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(B,2), median(B,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    %     plot(1:size(B,2), mode_B,'*','MarkerEdgeColor',c(1,:));
    ax = gca;
    ax.XTick = 1:n;
    ax.XTickLabel = xl;
    ax.XLim = [0.5 n+0.5];
    ax.YLabel.String = 'min';
    ax.Title.String = sprintf('B (mean=%G min)',round(mean_B(1),3,'significant'));
    ax.FontSize = 8;
    %      ax.YScale = 'log';
    %
    subplot(2,5,8), hold on,
    v = exp([out.v out.v_rep]);
    %     mode_v = mode(round(v,3,'significant'));
    mean_v = mean(v,1);
    plot(1:size(v,2), mean_v,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    plot(1:size(v,2), quantile(v,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(v,2), quantile(v,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(v,2), median(v,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    %     plot(1:size(v,2), mode_v,'*','MarkerEdgeColor',c(1,:));
    ax = gca;
    ax.XTick = 1:n;
    ax.XTickLabel = xl;
    ax.XLim = [0.5 n+0.5];
    %     ax.YLabel.String = 'a.u.';
    ax.Title.String = sprintf('v (mean=%G)',round(mean_v(1),3,'significant'));
    ax.FontSize = 8;
    %     ax.YScale = 'log';
    %
    %     subplot(2,5,10), hold on,
    %     E = [out.E out.E_rep];
    %     mode_E = mode(round(E,3,'significant'));
    %     mean_E = mean(E,1);
    %     plot(1:size(E,2), mean_E,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    %     plot(1:size(E,2), quantile(E,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    %     plot(1:size(E,2), quantile(E,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    %     plot(1:size(E,2), median(E,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    % %     plot(1:size(v,2), mode_v,'*','MarkerEdgeColor',c(1,:));
    %     ax = gca;
    %     ax.XTick = 1:n;
    %     ax.XTickLabel = xl;
    %     ax.XLim = [0.5 n+0.5];
    % %     ax.YLabel.String = 'a.u.';
    %     ax.Title.String = sprintf('E (mean=%0.2e min^{-1})',round(mean_E(1),3,'significant'));
    %     ax.FontSize = 8;
    subplot(2,5,9), hold on,
    Q = [out.Q out.Q_rep];
    M = Q + B.*log2(v);
    %     mode_M = mode(round(M,3,'significant'));
    mean_M = mean(M,1);
    plot(1:size(M,2), mean_M,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    plot(1:size(M,2), quantile(M,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(M,2), quantile(M,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(M,2), median(M,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    %     plot(1:size(M,2), mode_M,'*','MarkerEdgeColor',c(1,:));
    ax = gca;
    ax.XTick = 1:n;
    ax.XTickLabel = xl;
    ax.XLim = [0.5 n+0.5];
    ax.YLabel.String = 'min';
    ax.Title.String = sprintf('M (mean=%G min)',round(mean_M(1),3,'significant'));
    ax.FontSize = 8;
    %      ax.YScale = 'log';
    %
    subplot(2,5,10), hold on,
    D = [out.min_relative_doubling out.min_relative_doubling_rep];
    mean_D = mean(D,1);
    plot(1:size(D,2), mean_D,'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    plot(1:size(D,2), quantile(D,0.1,1),'^','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(D,2), quantile(D,0.9,1),'v','MarkerEdgeColor','none','MarkerFaceColor',c(2,:));
    plot(1:size(D,2), median(D,1),'d','MarkerEdgeColor','none','MarkerFaceColor',c(3,:));
    %     plot(1:size(M,2), mode_M,'*','MarkerEdgeColor',c(1,:));
    ax = gca;
    ax.XTick = 1:n;
    ax.XTickLabel = xl;
    ax.XLim = [0.5 n+0.5];
    ax.YLabel.String = 'min';
    ax.Title.String = sprintf('D (mean=%G min)',round(mean_D(1),3,'significant'));
    ax.FontSize = 8;
    
    %     subplot(2,5,10), hold on,
    %     mean_temp = nanmean(dataOD.temperature(ind));
    %     plot(dataOD.time(ind), dataOD.temperature(ind),'o','MarkerEdgeColor','none','MarkerFaceColor',c(1,:));
    %     ax = gca;
    %     ax.YLabel.String = '^{\circ}C';
    %     ax.Title.String = sprintf('Temp (mean=%G%s)',round(mean_temp,3,'significant'),' ^{\circ}C');
    %     ax.FontSize = 8;
    %     ax.XLim = [min(dataOD.time(ind)) max(dataOD.time(ind))];
    %     xlabel('Time (min)');
    %% write summary file
    fitSummary = fitOD.print();
    ind = cellfun(@isempty,strfind(fitSummary,'_pred')) & ...
        cellfun(@isempty,strfind(fitSummary,'***'));
    fitSummary = fitSummary(ind);
    fid = fopen(strcat(baseName,'_',samples{k},'_summary.txt'),'wt');
    for i =1:numel(fitSummary)
        fprintf(fid,'%s\n',fitSummary{i});
    end
    fclose(fid);
    %     writetable(fitSummaryTable,strcat(baseName,'_',samples{k},'_summary.txt'));
    %% save plots
    set(h,'PaperPositionMode','auto');
    print(h,strcat(baseName,'_',samples{k},'_plots.png'),'-dpng','-r0','-noui');
    print(h,strcat(baseName,'_',samples{k},'_plots.svg'),'-dsvg','-r0','-noui');
    crop(strcat(baseName,'_',samples{k},'_plots.png'));
    %%
    delete('*-data.R',strcat('*-output-*.csv'));
end
end
