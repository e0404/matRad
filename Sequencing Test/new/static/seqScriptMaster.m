dir = pwd;

%{
cd('TG119')
seqScript
cd(dir)

cd('Prostate')
seqScript
cd(dir)
%}
cd('H&N')
seqScript2deg
cd(dir)