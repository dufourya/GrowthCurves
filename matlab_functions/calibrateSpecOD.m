function calibrateSpecOD(DilutionTable)

[~,name,~] = fileparts(DilutionTable);
data = readtable(DilutionTable);
data = table2array(data);
rep = data(:,1);
[~,~,rep] = unique(rep);
dilution = data(:,2);
od = data(:,3);

referenceOD = 0.2;
%%

opts = optimoptions('lsqnonlin','Display','off',...
    'MaxFunctionEvaluations',10000,'MaxIterations',10000);

% f = @(p) (((exp(p(2)*(od.^p(4)-p(3))/(p(1)+p(3)-od.^p(4)))-1)-(dilution))./od).^2;

% f = @(p) (((exp(p(2)*(od-p(3))/(p(1)+p(3)-od))-1)-(dilution))./od).^2;
%f = @(p) ((abs(((od-p(3))./(p(1)-od)).^(p(2)))-dilution)./od).^2;

f = @(p) (od-(p(3) + p(1) .* ((dilution .^ p(2)) ./ (dilution .^ p(2) + p(4)))));

% pinit = ones(max(rep)+4,1);
pinit(1) = 1.35;
pinit(2) = 0.666;
pinit(3) = 0.0831;
pinit(4) = 1;
% pinit(4) = 5;
% plow = zeros(max(rep)+4,1);
% plow(3) = -Inf;
% plow(4) = 1;
% phigh = Inf * ones(max(rep)+4,1);
% phigh(1) = 1;
% pout = lsqnonlin(f,pinit,plow,phigh,opts);
pout = lsqnonlin(f,pinit,[0 0 0 0],[Inf Inf Inf Inf],opts);
%pout = [1.35 0.666 0.0831];
%%

h = figure;
c = lines(max(rep));
% OD = @(x) (pout(1).*x + pout(2).*(x).^pout(4)+pout(3));
% OD = @(x) (exp(pout(2)*(x.^pout(4)-pout(3))./(pout(1)+pout(3)-x.^pout(4)))-1);
%OD = @(x) abs((x-pout(3))./(pout(1)-x)).^(pout(2));
% OD = @(x) (exp(pout(2)*(x-pout(3))./(pout(1)+pout(3)-x))-1);
OD = @(x) pout(3) + pout(1) .* (x .^ pout(2) ./ (x .^ pout(2) + pout(4)));

plot(min(dilution):0.01:max(dilution)*1.1,OD(min(dilution):0.01:1.1*max(dilution)),'linewidth',1,'color','k');
hold on;
% plot([0 max(od)*1.1],[pout(3)./OD(referenceOD).*referenceOD max(od)*1.1+pout(3)./OD(referenceOD).*referenceOD],':k');
scatter(dilution, od, 50,c(rep,:),'filled'); hold off;

ylim([0 1.1*(max(od))]);
xlim([0 max(dilution)*1.1]);
axis square
legend({'Model','Data'},'Location','NorthWest');
xlabel('Relative cell concentration');
ylabel('Measured OD');

%%
s = sprintf('\nModel equation:\n\tMeasured OD = blank + \nalpha * (cells^{beta}) / (cells^{beta} + offset)\n\nFitted parameter values:\n\talpha = %.5g (Maximum OD)\n\tbeta = %.5g (non linear factor)\n\tblank = %.5g (blank OD)\n\toffset = %.5g',pout);
text(0.4*max(dilution),max(od)*0.25,s);
print(h,strcat(name,'.png'),'-dpng','-r300','-noui');

end