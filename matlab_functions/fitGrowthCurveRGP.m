function fitResults = fitGrowthCurveRGP(filename, varargin)

[~,basename,ext] = fileparts(filename);

%%
p = inputParser;
addRequired(p,'filename',@ischar);
addParameter(p,'samples','',@ischar);
addParameter(p,'CI',0.95,@isnumeric);
addParameter(p,'Save','');
addParameter(p,'Exclude','');
addParameter(p,'Calibration','none');%@(x) ismember(x,{'ecoli','vibrio'})
addParameter(p,'SkipPrint','');
addParameter(p,'DropTimePoints',0,@isnumeric);

parse(p,filename,varargin{:});
samples = p.Results.samples;
conf = p.Results.CI;
saveFormat = p.Results.Save;
excluded = p.Results.Exclude;
calibration = p.Results.Calibration;
skipPrint = p.Results.SkipPrint;
dropPoint = p.Results.DropTimePoints;
%%
switch lower(ext)
    case '.asc'
        % read data from ASC file (CSV)
        data = readtable(p.Results.filename, 'Format', 'auto', 'Filetype','text', 'ReadVariableNames',true,...
            'ReadRowNames',true,'MultipleDelimsAsOne',false);
        ind = find(strcmp(data.Properties.VariableNames,'RawData'));
        samplingtime = data(1,ind:end-1);
        samplingtime = str2double(strrep(samplingtime.Variables,'s',''))/60;
        temperature = data(2,ind:end-1);
        temperature = str2double(regexprep(temperature.Variables,'[^\d\.]',''));
        rawdata = data(3:end,ind:end-1);
        rawdata = cellfun(@str2double,rawdata.Variables);
        if ~strcmp('none',calibration)
            [caldata, alpha, beta, offset] = calibrateOD(rawdata, calibration);
        end
        metadata = data(3:end,1:ind-1);
        sample = string(repmat(metadata.Layout,1,size(caldata,2)));
        position = string(repmat(metadata.Properties.RowNames,1,size(caldata,2)));
        samplingtime = repmat(samplingtime,size(caldata,1),1);
        temperature = repmat(temperature,size(caldata,1),1);
        data = struct('sample',sample(:),...%'replicate',replicate(:),...
            'time',samplingtime(:),'OD',caldata(:),'temperature',...
            temperature(:),'position', position(:), 'ODraw', rawdata(:));
        data = struct2table(data);
    case '.txt'
        % read data from TXT file (TAB)
        data = readtable(p.Results.filename,'ReadVariableNames',true);
        data = table2struct(data,'ToScalar',true);
        data.sample = string(data.sample);
        if ~isfield(data,'position')
            data.position = string(data.replicate);
        else
            data.position = strcat(data.replicate,string(data.position));
        end
        if ~isfield(data,'ODraw')
            data.ODraw = data.OD;
        end
        data = struct2table(data);
        alpha = 1;
        beta = 1;
	offset = 1;
    otherwise
        error('File format not recognized');
end
%%

times = unique(data.time);

ind = data.time>=times(dropPoint+1);
data = data(ind,:);

datablank = data(contains(data.sample,'BL'),:);
ind = contains(data.sample,samples);
data = data(ind,:);
data = data(~ismember(data.position, excluded),:);
data.OD = data.OD - min(0,min(data.OD));

[wells,ia,~] = unique(data.position);
wells = wells(~ismember(wells,excluded));

newdata = [];
fitParams = [];
opts = optimoptions('lsqnonlin','Display','off');

fitResults = zeros(numel(wells),3);

fittingData = [];
fittingData = [fittingData; datablank];

figure,
for k = 1:numel(wells)
    %%
    ind = strcmp(data.position,wells{k});
    subdata = data(ind,:);
    [~,order] = sort(subdata.time);
    subdata = subdata(order,:);
    %subdata = subdata(order(2:end-1),:);
    
    %rgp = fitrgp(subdata.time,subdata.OD,'KernelFunction','squaredexponential','ConstantSigma',true,'Sigma',0.00001,'FitMethod','none');
    %ypred = resubPredict(rgp);
    %ypred = expconv2(subdata.time,subdata.OD,20);
    
    dOD = diff(subdata.OD)./diff(subdata.time);
    dOD = ([NaN; dOD] + [dOD; NaN])/2;
    
    ddOD = diff(dOD)./diff(subdata.time);
    ddOD = ([NaN; ddOD] + [ddOD; NaN])/2;
    
    dODraw = diff(subdata.ODraw)./diff(subdata.time);
    dODraw = ([NaN; dODraw] + [dODraw; NaN])/2;
    
    ddODraw = diff(dODraw)./diff(subdata.time);
    ddODraw = ([NaN; ddODraw] + [ddODraw; NaN])/2;
    
    idn = find(dODraw >=0 & ddODraw >= 0);
    idn = idn(find(diff(idn)==1,1,'first'));
    
    %     ldOD = dOD./subdata.OD;
    %         idn = find(dOD >=0 & ddOD >= 0, 1, );
    idx = find(ddOD(idn:end)<0 & dODraw(idn:end)>(0.5*max(dODraw(idn:end))),1,'first');
    idx = idx + idn - 2;
    r = dOD(idx);
    
    
%     idn = 1;
%     k
%     figure,
%     subplot(6,1,1), plot(subdata.ODraw); hold on; scatter([idn idx], [subdata.ODraw(idn) subdata.ODraw(idx)]);
%     subplot(6,1,2), plot(dODraw); hold on; scatter([idn idx], [dODraw(idn) dODraw(idx)]);
%     subplot(6,1,3), plot(ddODraw); hold on; scatter([idn idx], [ddODraw(idn) ddODraw(idx)]);
%     subplot(6,1,4), plot(subdata.OD); hold on; scatter([idn idx], [subdata.OD(idn) subdata.OD(idx)]);
%     subplot(6,1,5), plot(dOD); hold on; scatter([idn idx], [dOD(idn) dOD(idx)]);
%     xlim([0 numel(dOD)]);
%     subplot(6,1,6), plot(ddOD); hold on; scatter([idn idx], [ddOD(idn) ddOD(idx)]);
%     xlim([0 numel(dOD)]);
%     pause;
    
    if (idx-idn)>1
        f = @(p) ((p(1) + alpha * (p(2) .* exp(p(3) .* subdata.time(idn:idx))).^beta ./...
            ((p(2) .* exp(p(3) .* subdata.time(idn:idx)).^beta) + offset)) - (subdata.ODraw(idn:idx))) ./ subdata.ODraw(idn:idx);
        pout = lsqnonlin(f,[subdata.ODraw(idn) subdata.ODraw(idn)/10 2*r],[0 0 0],[Inf Inf Inf],opts);
    else
        error("Failed to fit well: %s. Check your data. This often occurs because dOD/OD has a large amount of noise early in the run.", wells{k});
    end
    fittingData = [fittingData; subdata(idn:idx, :)];
    temp = subdata(1,:);
    temp.blank = pout(1);
    temp.L0 = pout(2);
    temp.rate = pout(3);
    
    subdata.OD = ((subdata.ODraw - pout(1)) ./ (alpha + pout(1) - subdata.ODraw)).^ (1/beta);
    ii = find(imag(subdata.OD)~=0,1,'last');
    subdata.OD(1:ii) = NaN;
    subdata.OD = real(subdata.OD);
    subdata.log2OD = log2(subdata.OD);
    
    % rate
    dOD = diff(subdata.OD)./diff(subdata.time);
    subdata.dOD = ([dOD(1); dOD] + [dOD; dOD(end)])/2 ./subdata.OD;
    ii = find(subdata.dOD>=0,1,'first');
    subdata.OD(1:ii) = NaN;
    subdata.dOD(1:ii) = NaN;
    
    newdata = [newdata; subdata];
    fitParams = [fitParams; temp];
    
    fitResults(k,:) = pout;
    
    %                     ax = figure;
    
    scatter(subdata.time(idn:idx), (subdata.ODraw(idn:idx)));hold on;
    plot(subdata.time(idn:idx), (pout(1)+alpha*(pout(2).*exp(pout(3).*subdata.time(idn:idx))).^beta ./...
        ((pout(2).*exp(pout(3).*subdata.time(idn:idx)).^beta)+ offset))); %%hold off;
    ax.Children.YScale = 'log';
    xlabel('Time (min)');
    ylabel('Raw OD');
    hold off;
    saveas(gcf, sprintf('well_%s.png', wells{k}));
     clf
    %             pause;
    
end
fitResults = array2table(fitResults,'RowNames',{wells{:}},'VariableNames',{'Blank','InitialOD','GrowthRate'});
fitResults = addvars(fitResults,data.sample(ia),'NewVariableNames','Sample','Before','Blank');
if ismember('txt',saveFormat)
    writetable(fitResults,strcat(basename,'.summary.txt'),'WriteRowNames',true);
    writetable(fittingData,strcat(basename,'.fittingData.csv'),'WriteRowNames',true);
end

if (skipPrint == true)
    return;
end

%%
samples = unique(newdata.sample);
nb = regexp(samples,'[0-9]+','match');
[~,nb] = sort(str2double(cellfun(@(x) strjoin(x,''), nb, 'UniformOutput', false)));
samples = samples(nb);
%%
pd = makedist('Normal');
confint = icdf(pd,[(1-conf)/2 (1+conf)/2]);

summarydata = cell(numel(samples),3,2);

for k = 1:numel(samples)
    
    ind = strcmp(newdata.sample,samples{k}) & ~isnan(newdata.OD);
    params = fitParams(strcmp(fitParams.sample,samples{k}),:);
    data = newdata(ind,:);
    
    if size(data,1)>2
        
        rgpOD = fitrgp(data.time,data.OD,'KernelFunction','exponential','ConstantSigma',true,'Sigma',0.001,'FitMethod','none');
        predOD = resubPredict(rgpOD);
        
        % predOD = expconv2(data.time,data.OD,20);
        resOD = (data.OD-predOD).^2;
        [predresOD,neff1] = expconv2(data.time,resOD,60);
        [pos,~,pc] = unique(data.position);
        c = lines(numel(pos));
        
        h1 = figure('name',strcat(basename,'_',samples{k}),'position',[200 200 800 800]);
        
        [t, ii] = sort(data.time);
        
        ax1 = subplot(2,2,1); hold on;
        patch(ax1,[t; flipud(t)], [predOD(ii) + confint(1) * sqrt(predresOD(ii)./neff1); flipud(predOD(ii) + confint(2) * sqrt(predresOD(ii)./neff1))],[0.9 0.9 0.9],...
            'EdgeColor','none');
        r =[];
        for j = 1:numel(pos)
            r(j) = scatter(ax1,data.time(pc==j),data.OD(pc==j),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
        end
        plot(ax1,t,predOD(ii),'Linewidth',2,'Color',[0.3 0.3 0.3]);
        %         ylim([0 max(predOD(ii) + confint(2) * sqrt(predresOD(ii)))]);
        axis tight;
        xlim([0 max(t)]);
        axis square;
        xlabel('Time (min)');
        ylabel('OD_{600}');
        title(ax1,strcat(basename,'_',samples{k}),'Interpreter','none');
        legend(ax1,r,cellstr(pos),'Location','NorthWest','Fontsize',10);
        ax1.FontSize = 12;
        summarydata{k,1,1} = [[t; flipud(t)], [predOD(ii) + confint(1) * sqrt(predresOD(ii)./neff1(ii)); flipud(predOD(ii) + confint(2) * sqrt(predresOD(ii)./neff1(ii)))]];
        summarydata{k,1,2} = [t,predOD(ii)];
        
        indx = ~isnan(data.dOD) & data.time>max([0; data.time(data.OD<0.1*max(data.OD) & (data.dOD>4*mean(params.rate) | data.dOD<0.25*mean(params.rate)))]);
        
        if sum(indx)>1
            preddOD = nan(size(data,1),1);
            rgpdOD = fitrgp(data.time(indx),data.dOD(indx),'KernelFunction','squaredexponential','ConstantSigma',true,'Sigma',0.0001,'FitMethod','none');
            preddOD(indx) = resubPredict(rgpdOD);
            
            %             preddOD(indx) = expconv2(data.time(indx),data.dOD(indx),20);
            
            resdOD = (data.dOD-preddOD).^2;
            [predresdOD,neff2] = expconv2(data.time,resdOD,60);
            
            [~,imaxOD] = max(predOD);
            
            indg = preddOD>0 & data.time<=data.time(imaxOD) & data.time>max([0; data.time(data.OD<0.1*max(data.OD) & (data.dOD>4*mean(params.rate) | data.dOD<0.25*mean(params.rate)))]);
            
            if sum(indg)>0
                pODdOD = nan(size(data,1),1);
                rgpODdOD = fitrgp(log2(data.OD(indg)),data.dOD(indg),'KernelFunction','squaredexponential','ConstantSigma',true,'Sigma',0.0001,'FitMethod','none');
                pODdOD(indg) = resubPredict(rgpODdOD);
                % pODdOD(indg) = expconv2(log2(data.OD(indg)),data.dOD(indg),0.2);
                resODdOD = (data.dOD-pODdOD).^2;
                [predresODdOD,neff3] = expconv2(log2(data.OD),resODdOD,0.25);
                
                ax2 = subplot(2,2,3); hold on;
                iii = ii(~isnan(preddOD(ii)));
                hp = patch(ax2,[data.time(iii); flipud(data.time(iii))], [preddOD(iii) + confint(1) * sqrt(predresdOD(iii)./neff2(iii)); flipud(preddOD(iii) + confint(2) * sqrt(predresdOD(iii)./neff2(iii)))],[0.9 0.9 0.9],...
                    'EdgeColor','none');
                for j = 1:numel(pos)
                    r(j) = scatter(ax2,data.time(pc==j),data.dOD(pc==j),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
                    plot(ax2,[-1 max(t)+1], [1 1] * params.rate(strcmp(params.position,pos{j})),':','Color',c(j,:),'LineWidth',1);
                end
                %             scatter(ax2,data.time,data.dOD,20,categorical(data.position),'filled','MarkerFaceAlpha',0.5);
                plot(ax2,data.time(iii),preddOD(iii),'Linewidth',2,'Color',[0.3 0.3 0.3]);
                plot(ax2,[-1 max(t)+1], [0 0],'k');
                ax2.XLim = ax1.XLim;
                ax2.YLim = [min(preddOD(iii) + confint(1) * sqrt(predresdOD(iii)./neff2(iii))) 1.5*mean(params.rate)]; %max(preddOD(iii) + confint(2) * sqrt(predresdOD(iii)./neff2(iii)))];
                axis square;
                xlabel('Time (min)');
                ylabel('Growth rate (min^{-1})');
                ax2.FontSize = 12;
                legend(hp,{sprintf('%g conf. SE',conf)},'Location','Best','Interpreter','none','Fontsize',10);
                title(ax2,sprintf('OD calibration: %s',calibration));
                summarydata{k,2,1} = [[data.time(iii); flipud(data.time(iii))], [preddOD(iii) + confint(1) * sqrt(predresdOD(iii)./neff2(iii)); flipud(preddOD(iii) + confint(2) * sqrt(predresdOD(iii)./neff2(iii)))]];
                summarydata{k,2,2} = [data.time(iii),preddOD(iii)];
                
                [~, ii] = sort(data.OD);
                ax3 = subplot(2,2,4); hold on;
                iii = ii(preddOD(ii)>0 & predOD(ii)>0 & ~isnan(pODdOD(ii)));
                if sum(iii)>1
                    patch(ax3,log2([data.OD(iii); flipud(data.OD(iii))]), [pODdOD(iii) + confint(1) * sqrt(predresODdOD(iii)./neff3(iii)); flipud(pODdOD(iii) + confint(2) * sqrt(predresODdOD(iii)./neff3(iii)))],[0.9 0.9 0.9],...
                        'EdgeColor','none');
                    for j = 1:numel(pos)
                        r(j) = scatter(ax3,real(log2(data.OD(pc==j))),data.dOD(pc==j),20,'filled','MarkerFaceAlpha',0.6,'MarkerFaceColor',c(j,:));
                        plot(ax3,[-13 10], [1 1] * params.rate(strcmp(params.position,pos{j})),':','Color',c(j,:),'LineWidth',1);
                        
                    end
                    %                 scatter(ax3,real(log2(data.OD)),data.dOD,20,categorical(data.position),'filled','MarkerFaceAlpha',0.5);
                    plot(ax3,real(log2(data.OD(iii))),pODdOD(iii),'Linewidth',2,'Color',[0.3 0.3 0.3]);
                    plot(ax3,[-13 10], [0 0],'k');
                    %                     xlim([-7 max(log2(data.OD(iii)))]);
                    axis square;
                    ax3.YLim = [min(pODdOD(iii) + confint(1) * sqrt(predresODdOD(iii)./neff3(iii))) max(pODdOD(iii) + confint(2) * sqrt(predresODdOD(iii)./neff3(iii)))];
                    xlabel('OD_{600}');
                    ylabel('Growth rate (min^{-1})');
                    ax3.FontSize = 12;
                    ax3.XTick = -12:3;
                    ax3.XTickLabel = {'2^{-12}' '2^{-11}' '2^{-10}' '2^{-9}' '2^{-8}' '2^{-7}' '2^{-6}' '2^{-5}' '2^{-4}' '2^{-3}' '2^{-2}' '2^{-1}' '2^0' '2^1' '2^2' '2^3'};
                    ax3.XLim = [max(min(real(log2(data.OD(iii)))),-13) max(real(log2(data.OD(iii))))];
                    ax3.YLim(1) = min(ax3.YLim(1),ax2.YLim(1));
                    ax2.YLim(1) = ax3.YLim(1);
                    ax3.YLim(2) = 1.5*mean(params.rate); % max(ax3.YLim(2),ax2.YLim(2));
                    %                     ax2.YLim(2) = ax3.YLim(2);
                    
                    summarydata{k,3,1} = [log2([data.OD(iii); flipud(data.OD(iii))]), [pODdOD(iii) + confint(1) * sqrt(predresODdOD(iii)./neff3(iii)); flipud(pODdOD(iii) + confint(2) * sqrt(predresODdOD(iii)./neff3(iii)))]];
                    summarydata{k,3,2} = [log2(data.OD(iii)),pODdOD(iii)];
                end
            end
        end
        
    end
    
    %%
    
    ind = strcmp(fitParams.sample,samples{k}) & fitParams.rate > 0;
    ax4 = subplot(2,4,3); hold on;
    scatter(2*ones(sum(ind),1),fitParams.blank(ind),100,c,'filled','MarkerFaceAlpha',0.6);
    scatter(ones(numel(datablank.ODraw),1), datablank.ODraw,100,'k','filled','MarkerFaceAlpha',0.05);
    ax4.XLim = [0.5 2.5];
    ax4.XTick = [1 2];
    ax4.XTickLabel = {'Blanks','Wells'};
    ax4.YLabel.String = 'OD_{600}';
    
    title('Fitted blanks');
    ax5 = subplot(2,4,4);
    scatter(fitParams.L0(ind),fitParams.rate(ind),100,c,'filled','MarkerFaceAlpha',0.6);
    ax5.XLabel.String = 'Initial OD_{600}';
    ax5.YLabel.String = 'Initial growth rate';
    title('Fitted growth');
    
    %     yl = ax5.YLim;
    nticks = 15; %numel(ax5.YTick)+2;
    yl = [min(fitParams.rate(ind))*0.95, max(fitParams.rate(ind))*1.05];
    ax5.YLim = yl;
    ticks = fliplr(unique(floor(exp(linspace(log(log(2)/yl(1)),log(log(2)./yl(2)),nticks)))));
    
    yyaxis(ax5,'right');
    ax5.YLim = yl;
    ax5.YTick = log(2)./ticks;
    ax5.YTickLabel = num2str(ticks');
    ax5.YLabel.String = 'Doubling time (min)';
    
    %%
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
c = linspecer(numel(samples));
h2 = figure('name',basename,'position',[200 200 800 800]);
miny = Inf;
maxy = -Inf;
maxx = -Inf;
maxt = -Inf;
minx = Inf;

ax2 = subplot(2,2,3); hold on;
ax3 = subplot(2,2,4); hold on;
ax4 = subplot(2,4,3); hold on;
scatter(ax4, ones(numel(datablank.ODraw),1), datablank.ODraw,100,'k','filled','MarkerFaceAlpha',0.05);
ax5 = subplot(2,4,4); hold on;

for k = 1:numel(samples)
    ax1 = subplot(2,2,1); hold on;
    patch(ax1, summarydata{k,1,1}(:,1),summarydata{k,1,1}(:,2),c(k,:),'FaceAlpha',0.3,'EdgeColor','none');
    h(k) = plot(ax1, summarydata{k,1,2}(:,1),summarydata{k,1,2}(:,2),'Color',c(k,:),'LineWidth',2);
    %     idx = summarydata{k,1,2}(:,2)>2^-7;
    if ~isempty(summarydata{k,3,1})
        patch(ax2, summarydata{k,2,1}(:,1),summarydata{k,2,1}(:,2),c(k,:),'FaceAlpha',0.3,'EdgeColor','none');
        plot(ax2, summarydata{k,2,2}(:,1),summarydata{k,2,2}(:,2),'Color',c(k,:),'LineWidth',2);
        plot(ax2, [-1 10000],[1 1] * mean(fitParams.rate(strcmp(fitParams.sample,samples{k}))),':','Color',c(k,:),'LineWidth',1);
        title(ax2,sprintf('OD calibration: %s',calibration));
        
        patch(ax3, summarydata{k,3,1}(:,1),summarydata{k,3,1}(:,2),c(k,:),'FaceAlpha',0.3,'EdgeColor','none');
        plot(ax3, summarydata{k,3,2}(:,1),summarydata{k,3,2}(:,2),'Color',c(k,:),'LineWidth',2);
        plot(ax3, [-13 10],[1 1] * mean(fitParams.rate(strcmp(fitParams.sample,samples{k}))),':','Color',c(k,:),'LineWidth',1);
        miny = min([miny; summarydata{k,2,1}(:,2); summarydata{k,3,1}(:,2)]);
        maxy = max(maxy, mean(fitParams.rate(strcmp(fitParams.sample,samples{k}))));
        minx = min([minx; summarydata{k,3,2}(:,1)]);
        maxx = max([maxx; summarydata{k,3,2}(:,1)]);
        maxt = max([maxt; summarydata{k,1,1}(:,1)]);
        
        ind = strcmp(fitParams.sample,samples{k}) & fitParams.rate > 0;
        scatter(ax4,2*ones(sum(ind),1),fitParams.blank(ind),100,c(k,:),'filled','MarkerFaceAlpha',0.6);
        scatter(ax5,fitParams.L0(ind),fitParams.rate(ind),100,c(k,:),'filled','MarkerFaceAlpha',0.6);
        
    end
end


title(ax1,basename,'Interpreter','none');
ax1.XLabel.String = 'Time (min)';
ax1.YLabel.String = 'OD_{600}';
ax1.FontSize = 12;
ax1.PlotBoxAspectRatio = [1 1 1];
ax1.XLim = [0 maxt];
ax1.YLim(1) = 0;
legend(h,samples,'Location','NorthWest','Interpreter','none','Fontsize',10);

ax2.XLabel.String = 'Time (min)';
ax2.YLabel.String = 'Growth rate (min^{-1})';
ax2.FontSize = 12;
ax2.YLim = [max(miny,-0.5*maxy) 1.5*maxy];
ax2.XLim = ax1.XLim;
ax2.PlotBoxAspectRatio = [1 1 1];
plot(ax2,ax2.XLim, [0 0],'k');
hp = patch(ax2,0,0,[0.9 0.9 0.9],'EdgeColor','none');
legend(hp,{sprintf('%g conf. SE',conf)},'Location','Best','Interpreter','none','Fontsize',10);
% axes(ax2,'YAxisLocation','right');

% ax2.YLim(2) = max(ax2.YLim(2), ax3.YLim(2));
ax3.XLabel.String = 'OD_{600}';
ax3.YLabel.String = 'Growth rate (min^{-1})';
ax3.FontSize = 12;
ax3.FontSize = 12; ax3.XTick = -12:3;
ax3.XTickLabel = {'2^{-12}' '2^{-11}' '2^{-10}' '2^{-9}' '2^{-8}' '2^{-7}' '2^{-6}' '2^{-5}' '2^{-4}' '2^{-3}' '2^{-2}' '2^{-1}' '2^0' '2^1' '2^2' '2^3'};
ax3.XLim(1) = max(-13,minx);
ax3.YLim = ax2.YLim;
ax3.XLim(2) = maxx;
ax3.PlotBoxAspectRatio = [1 1 1];
plot(ax3,ax3.XLim, [0 0],'k');

ax5.XLabel.String = 'Initial OD_{600}';
ax5.YLabel.String = 'Initial growth rate (min^{-1})';
%     axis square;
ax5.Title.String = 'Fitted growth';
% ax5.PlotBoxAspectRatio = [1 2 1];
yl = ax5.YLim;
nticks = 20; %numel(ax5.YTick)+2;

ax4.XLim = [0.5 2.5];
ax4.XTick = [1 2];
ax4.XTickLabel = {'Blanks','Wells'};
ax4.YLabel.String = 'OD_{600}';
%     axis square;
ax4.Title.String = 'Fitted blanks';
% ax4.PlotBoxAspectRatio = [1 2 1];

yyaxis(ax5,'right');
ax5.YLim = yl;
yl(1) = max(yl(1), 0.0025);
ticks = fliplr(unique(floor(exp(linspace(log(log(2)/yl(1)),log(log(2)./yl(2)),nticks)))));
ax5.YTick = log(2)./ticks;
ax5.YTickLabel = num2str(ticks');
ax5.YLabel.String = 'Doubling time (min)';


if ismember('svg',saveFormat)
    set(h2,'PaperPositionMode','auto');
    print(h2,sprintf('%s_Summary.svg',basename),'-dsvg','-r300','-noui');
end
if ismember('png',saveFormat)
    set(h2,'PaperPositionMode','auto');
    print(h2,sprintf('%s_Summary.png',basename),'-dpng','-r300','-noui');
    crop(sprintf('%s_Summary.png',basename));
end

end

function [out, neff] = expconv2(t,data,lag)
out = nan(numel(t),1);
neff = nan(numel(t),1);

ind = ~isnan(data);
t = t(ind);
data = data(ind);
temp = nan(numel(t),1);
tempn = nan(numel(t),1);

for k = 1:numel(t)
    tt = t-t(k);
    w = exp(-tt.^2/(2*lag^2));
    tempn(k) = sum(w);
    w = w/sum(w);
    temp(k) = sum(data.*w);
end

out(ind) = temp;
neff(ind) = tempn;
end

function [calibrated, a, b, k] = calibrateOD(data,calibration)
switch calibration
    case 'salmonellaVB2YE'
        %TH9677 in VB+2%YE (w/v) from Josh's 20190215 stationary phase cal
        a = 2.6451;
        b = 1.0711;
        c = 0.108;
        k = 0.98087;
        calibrated = (((data .* k) - (c * k)) ./ (a + c - data)).^ (1/b);
        calibrated(imag(calibrated) ~= 0) = nan;
        calibrated = real(calibrated);
        
    case 'classicalM9Pyr'
        %Classical strain in M9Pyruvate from Nhu's 20190206 stationary phase cal
        a = 24.024;
        b = 0.99579;
        c = 0.067021;
        k = 42.784;
        calibrated = (((data .* k) - (c * k)) ./ (a + c - data)).^ (1/b);
        calibrated(imag(calibrated) ~= 0) = nan;
        calibrated = real(calibrated);
        
    case 'LTEE'
        %Classical strain in M9Pyruvate from Nhu's 20190206 stationary phase cal
        a = 24.024;
        b = 0.99579;
        c = 0.08;
        k = 42.784;
        calibrated = (((data .* k) - (c * k)) ./ (a + c - data)).^ (1/b);
        calibrated(imag(calibrated) ~= 0) = nan;
        calibrated = real(calibrated);
        
    otherwise
        error('No calibration available');
end
end
