function runDEMO1
% delete previous .txt files:
layer_name = {'joint'};
for i = 1:length(layer_name)
    delete([layer_name{i} '*.txt']);
end
system('PDSTEP_demo.exe');
graphJoint