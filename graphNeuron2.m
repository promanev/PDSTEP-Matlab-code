function graphNeuron2
% works with a single file
system('PDSTEP_demo.exe');
layer_name = {'Input', 'Hidden', 'Output'};
layer_size=[1,2;3,4;5,16];

fID = fopen('neuron.txt','r');
a = fscanf(fID,'%f',[16 Inf])';
fclose(fID);  

% load touches:
fID = fopen('swingFootTouch.txt','r');
b = fscanf(fID,'%f',[1 Inf])';
fclose(fID);
xVals = 1:length(b);
xVals = xVals*10;

mainfig = figure('units','normalized','outerposition',[0 0 1 1]);
for i = 1:length(layer_name)
    subplot(1,length(layer_name),i)
    
    for n=layer_size(i,1):layer_size(i,2)
        
        plot(a(:,n))
        hold all
        grid on
        if n<3
            count = n;
        elseif and(n>=3,n<5)
            count = n-2;
        else
            count = n-4;
        end
        %disp(['Using ' num2str(n) '-th column to plot ' layer_name{i} ' #' num2str(count)])
        legend_text{count} = [layer_name{i} ' neuron ' int2str(count)];
    end
    
    ylimits = ylim;
    line([length(a(:,n))/2 length(a(:,n))/2],[ylimits(1) ylimits(2)])
    
    if i==1
        figure
        plot(a(:,1),'-o')
        hold on
        plot(xVals,b,'-x')
        hold off
        figure(mainfig)
    end
    legend(legend_text{:},'Location', 'Best')
    clear('legend_text');
    title([layer_name{i} ' Layer'])
end