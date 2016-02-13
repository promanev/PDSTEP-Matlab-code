function runDEMO
% delete previous .txt files:
layer_name = {'Input', 'Hidden', 'Output'};
for i = 1:length(layer_name)
    delete([layer_name{i} '*.txt']);
end
system('PDSTEP_demo.exe')
graphNeuron