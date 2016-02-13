function graphNeuron
layer_name = {'Input', 'Hidden', 'Output'};
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(layer_name)
    input_list = dir([layer_name{i} '*.txt']);
    input_list_cell = struct2cell(input_list);
    sz = size(input_list_cell); 

    subplot(1,length(layer_name),i)
    for n=1:sz(2)
        fID = fopen([layer_name{i} int2str(n) '.txt'],'r');
        a = fscanf(fID,'%f',[2 Inf]);
        plot(a(1,:),a(2,:))
        hold on
        fclose(fID);
        legend_text{n} = [layer_name{i} ' neuron ' int2str(n)];
    end
    legend(legend_text{:},'Location', 'Best')
    title([layer_name{i} ' Layer'])
end