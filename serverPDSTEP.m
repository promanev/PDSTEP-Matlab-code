
function serverPDSTEP(output_port, number_of_retries)
% SERVER Write a message over the specified port
% 
% Usage - server(message, output_port, number_of_retries)
import java.net.ServerSocket
import java.io.*

if nargin < 2
    number_of_retries = 20; % set to -1 for infinite
end

if nargin < 1
    output_port = 27015; % port hardcoded in C++ code
end

retry             = 0;

server_socket  = [];
output_socket  = [];

weights = [0.1 0.2222 0.3333333 44444.44 55.5555];
% while true
% 
%     retry = retry + 1;

try
%         if ((number_of_retries > 0) && (retry > number_of_retries))
%             fprintf(1, 'Too many retries\n');
%             break;
%         end

%     fprintf(1, ['Try %d waiting for client to connect to this ' ...
%                 'host on port : %d\n'], retry, output_port);



    % wait for 1 second for client to connect server socket
%     disp('Above server socket thing')
    server_socket = ServerSocket(output_port);
    
%     start = tic;
    server_socket.setSoTimeout(500);
%     disp(['Elapsed time ', num2str(toc(start))])
    
    system('PDSTEP_demo.exe &');
%     disp('I am below .exe')
    input_socket = server_socket.accept;
%     buffer_size = server_socket.getReceiveBufferSize;
%     disp(['Buffer size is: ', int2str(buffer_size)])
    input_stream = input_socket.getInputStream;
    d_input_stream = DataInputStream(input_stream);
% Hardcoded the length of fitness value in bytes = 8:
    message = zeros(1, 8, 'uint8');
    for i = 1:8
        message(i) = d_input_stream.read;
    end

    message = char(message);
    disp(['Recevied this value: ' message])
% OUTPUT PART:       
    output_stream   = input_socket.getOutputStream;
    d_output_stream = DataOutputStream(output_stream);
% 
    fprintf(1, 'READY TO OUTPUT\n');
%     pause(0.5)
%     for i=1:length(weights)
%         % output the data over the DataOutputStream
%         % Convert to stream of bytes
%         sendMessage = num2str(weights(i));
%         disp(['Inside writing cycle, step ', int2str(i)])
%         fprintf(1, 'Writing %d bytes\n', length(sendMessage));
%         disp(['Sending this value ' sendMessage])
%         writer = OutputStreamWriter(output_stream);
%         buffwriter = BufferedWriter(writer);
%         buffwriter.write(sendMessage,0,length(sendMessage));
%         buffwriter.flush;
% %         d_output_stream.writeBytes(sendMessage);
% %         d_output_stream.flush;
% %         pause(0.4)
%     end
    % Convert to stream of bytes
    sendMessage = num2str(weights);
    disp(['Sending this string stream: ' sendMessage])
    pause(1);
    d_output_stream.writeBytes(sendMessage);
    d_output_stream.flush;
%     writer.write(sendMessage,0,length(sendMessage));
%     writer.flush;
    pause(1);
    % clean up
    disp('Message sent, closing sockets')
    server_socket.close;
    input_socket.close; 
%         break;

catch
    if ~isempty(server_socket)
        server_socket.close
    end

    if ~isempty(output_socket)
        output_socket.close
    end

    % pause before retrying
%         pause(1);
end
% end
