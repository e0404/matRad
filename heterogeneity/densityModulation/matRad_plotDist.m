

density = 1.05;
numOfLungVoxels = 1e6;
n1 = 2.5; p1 = 0.25;
a = p1*(n1-1);
b = (1-p1)*(n1-1);
samples = betaincinv(rand([numOfLungVoxels,1]),a,b)*density;

n2 = 5; p2 = 0.6;
a = p2*(n2-1);
b = (1-p2)*(n2-1);
samples2 = betaincinv(rand([numOfLungVoxels,1]),a,b)*density;

n3 = 10; p3 = 0.9;
a = p3*(n3-1);
b = (1-p3)*(n3-1);
samples3 = betaincinv(rand([numOfLungVoxels,1]),a,b)*density;

% samples = sum(rand([numOfLungVoxels,round(n)]) < p, 2);
% samples = samples ./ round(n);

edges = linspace(0,1*density,41);
centers = (edges(1:end-1)+edges(2:end))/2;


f2 = figure;
h = histogram(samples,edges);
hold on
h2 = histogram(samples2,edges);
h3 = histogram(samples3,edges);
h.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
ylabel('probability')
xlabel('density [g/cm^3]')

legend({['n=',num2str(n1),' p=',num2str(p1)],['n=',num2str(n2),' p=',num2str(p2)],['n=',num2str(n3),' p=',num2str(p3)]},'Location','northeast')
% title(['Var = ',num2str(round(var(samples),4)),', mean = ',num2str(round(mean(samples),4))])
% matlab2tikz('pdfBeta.tex')


%%



density = 1.05;
numOfLungVoxels = 1e6;
n1 = 3; p1 = 0.25;
samples1 = sum(rand([numOfLungVoxels,round(n1)]) < p1, 2);
samples1 = samples1 ./ round(n1);

n2 = 5; p2 = 0.6;
samples2 = sum(rand([numOfLungVoxels,round(n2)]) < p2, 2);
samples2 = samples2 ./ round(n2);


n3 = 10; p3 = 0.9;
samples3 = sum(rand([numOfLungVoxels,round(n3)]) < p3, 2);
samples3 = samples3 ./ round(n3);


% samples = sum(rand([numOfLungVoxels,round(n)]) < p, 2);
% samples = samples ./ round(n);

edges = linspace(0,1*density,41);


f1 = figure;
h = histogram(samples1,edges);
hold on
h2 = histogram(samples2,edges);
h3 = histogram(samples3,edges);
h.Normalization = 'probability';
h2.Normalization = 'probability';
h3.Normalization = 'probability';
ylabel('probability')
xlabel('density [g/cm^3]')
ylim([0 0.5])

legend({['n=',num2str(n1),' p=',num2str(p1)],['n=',num2str(n2),' p=',num2str(p2)],['n=',num2str(n3),' p=',num2str(p3)]},'Location','northeast')
% title(['Var = ',num2str(round(var(samples),4)),', mean = ',num2str(round(mean(samples),4))])
matlab2tikz('pdfBino.tikz')