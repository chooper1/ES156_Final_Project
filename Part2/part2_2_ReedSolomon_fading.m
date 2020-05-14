%ES 156 Final Project
%Code for Part 2-2
%By: Coleman Hooper

%Parameters
%channel power constraint : P
%gaussian error standard deviation: sigma
%number of bits per transmission: num_bits
%number of simulation iterations: resolution
%scaling factor: rep_number

%E_b == P/num_bits for part 2-2

%this script was used for the simulations in part 2-2 with additive and
%multiplicative noise

clear all;

num_bits = 6; 
rep_number = 7/3; %used as a scaling parameter here 
resolution = 100;

%define R-S encoder and decoder
N = 7;  % Codeword length
K = 3;  % Message length
rsEncoder = comm.RSEncoder(N,K,'BitInput',false);
rsDecoder = comm.RSDecoder(N,K,'BitInput',false);

%testing different P and sigma values
P_lim = [0.1;0.1;10];
sig_lim = [0.1;0.1;10];
for P=P_lim(1):P_lim(2):P_lim(3) %varying P
    for sigma = sig_lim(1):sig_lim(2):sig_lim(3) %varying sigma
        %counters to keep track of total transmissions and errors
        total = 0;
        wrong = 0;
        for l = 1:resolution
            for k=1:num_bits % number of grouped transmissions 
                %generate random message
                for i=1:3
                    bits(i) = randi([0 1]);
                end
                input_bits(((k-1)*3+1):k*3) = bits(1:3);
                %encode the message to r-s code
                encData = rsEncoder(transpose(bits(1:3)));
                for i=1:7
                    coded(i) = mod(encData(i), 2);  
                end
                transmit(((k-1)*7+1):k*7) = coded;
            end
            num_coms = length(transmit) / num_bits; 
            for j=1:num_coms
                if num_bits > 1
                    %encode the r-s code to symbols
                    encoding(2*j-1:2*j) = encode_channel(transmit((j-1)*num_bits+1:j*num_bits), P, num_bits, rep_number);
                    %channel noise
                    n(1) = normrnd(0,sigma);
                    n(2) = normrnd(0,sigma);
                    h(1) = normrnd(0,1); 
                    h(2) = normrnd(0,1); 
                    encoding(2*j-1) = h(1)*encoding(2*j-1) +  n(1);
                    encoding(2*j) = h(2)*encoding(2*j) +  n(2);
                    %decoding the symbols to bits
                    bits_out_temp((j-1)*num_bits+1:j*num_bits) = decode_channel(encoding(2*j-1:2*j), P, num_bits, rep_number, h);
                else
                    %encode the r-s code to symbols
                    encoding(j) = encode_channel(transmit((j-1)*num_bits+1:j*num_bits), P, num_bits, rep_number);
                    %channel noise
                    n = normrnd(0,sigma);
                    h = normrnd(0,1);
                    encoding(j) = h*encoding(j) +  n;
                    %decoding the symbols to bits
                    bits_out_temp((j-1)*num_bits+1:j*num_bits) = decode_channel(encoding(j), P, num_bits, rep_number, h);
                end
            end
            for j = 1:(length(bits_out_temp)/7)
                %decode the output r-s message
                bits_out = rsDecoder(transpose(bits_out_temp(((j-1)*7+1):j*7)));
                total = total + 3;
                bit_check = input_bits(((j-1)*3+1):j*3);
                for i=1:3
                    if bit_check(i) ~= bits_out(i)
                        wrong = wrong + 1; 
                    end
                end
            end
        end
        %keeping track of BER and EPB over noise ratio for plotting
        P_e(round(P*10),round(sigma*10)) = log10(wrong/total);
        EPB_noise(round(P*10),round(sigma*10)) = 20*log10((P/num_bits)/(sigma^2));
    end
end

%plotting BER versus Noise Ratio (Simulated)
figure;
hold on;
xlabel("Energy-Per-Bit over Noise Ratio (dB)");
ylabel("Bit Error Probability");
for i=1:1:size(P_e,1)
    scatter(EPB_noise(i,:), P_e(i,:), 1, 'k');
end

%Encoding function
%this function scales the power to meet average power constraints,
%and then encodes the message into symbols
function codeword = encode_channel(bits, P, num_bits, rep_number)
    if num_bits == 1
        %scaling to meet average power constraints and encode the message as symbols
        if bits(1) == 1
            codeword(1) = sqrt(P/rep_number);
        else
            codeword(1) = -sqrt(P/rep_number);
        end
    elseif num_bits == 2
        %scaling to meet average power constraints and encode the message as symbols
        power_scaled = sqrt(P/(2*rep_number));
        if bits(1) == 1 
            codeword(1) = power_scaled;
        else
            codeword(1) = -power_scaled;
        end
        if bits(2) == 1
            codeword(2) = power_scaled;
        else
            codeword(2) = -power_scaled;
        end
    elseif num_bits == 4
        %scaling to meet average power constraints and encode the message as symbols
        power_scaled = sqrt(P/(10*rep_number));
        if bits(1) == 0 && bits(2) == 0
            codeword(1) = -3*power_scaled;
        elseif bits(1) == 0 && bits(2) == 1
            codeword(1) = -1*power_scaled;
        elseif bits(1) == 1 && bits(2) == 1
            codeword(1) = power_scaled;
        elseif bits(1) == 1 && bits(2) == 0
            codeword(1) = 3*power_scaled;
        end
        if bits(3) == 0 && bits(4) == 0
            codeword(2) = -3*power_scaled;
        elseif bits(3) == 0 && bits(4) == 1
            codeword(2) = -1*power_scaled;
        elseif bits(3) == 1 && bits(4) == 1
            codeword(2) = power_scaled;
        elseif bits(3) == 1 && bits(4) == 0
            codeword(2) = 3*power_scaled;
        end
    elseif num_bits == 6
        %scaling to meet average power constraints and encode the message as symbols
        power_scaled = sqrt(P/(42*rep_number));
        if bits(1) == 0 && bits(2) == 0 && bits(3) == 0 
            codeword(1) = -7*power_scaled;
        elseif bits(1) == 0 && bits(2) == 0 && bits(3) == 1
            codeword(1) = -5*power_scaled;
        elseif bits(1) == 0 && bits(2) == 1 && bits(3) == 1 
            codeword(1) = -3*power_scaled;
        elseif bits(1) == 0 && bits(2) == 1 && bits(3) == 0 
            codeword(1) = -1*power_scaled;
        elseif bits(1) == 1 && bits(2) == 1 && bits(3) == 0 
            codeword(1) = 1*power_scaled;
        elseif bits(1) == 1 && bits(2) == 1 && bits(3) == 1 
            codeword(1) = 3*power_scaled;
        elseif bits(1) == 1 && bits(2) == 0 && bits(3) == 1 
            codeword(1) = 5*power_scaled;
        elseif bits(1) == 1 && bits(2) == 0 && bits(3) == 0 
            codeword(1) = 7*power_scaled;
        end
        
        if bits(4) == 0 && bits(5) == 0 && bits(6) == 0 
            codeword(2) = -7*power_scaled;
        elseif bits(4) == 0 && bits(5) == 0 && bits(6) == 1
            codeword(2) = -5*power_scaled;
        elseif bits(4) == 0 && bits(5) == 1 && bits(6) == 1 
            codeword(2) = -3*power_scaled;
        elseif bits(4) == 0 && bits(5) == 1 && bits(6) == 0 
            codeword(2) = -1*power_scaled;
        elseif bits(4) == 1 && bits(5) == 1 && bits(6) == 0 
            codeword(2) = 1*power_scaled;
        elseif bits(4) == 1 && bits(5) == 1 && bits(6) == 1 
            codeword(2) = 3*power_scaled;
        elseif bits(4) == 1 && bits(5) == 0 && bits(6) == 1 
            codeword(2) = 5*power_scaled;
        elseif bits(4) == 1 && bits(5) == 0 && bits(6) == 0 
            codeword(2) = 7*power_scaled;
        end
    end
end

%Decoding function
%this function scales the power to meet average power constraints,
%and then decodes the symbols into bits
function bits = decode_channel(codeword, P, num_bits, rep_number, h)
    %accounting for multiplicative noise
    if num_bits == 1
        codeword = codeword / h; 
    else
        codeword(1) = codeword(1) / h(1);
        codeword(2) = codeword(2) / h(2);
    end
    
    if num_bits == 1
        %scaling to meet average power constraints and decoding the symbols
        %to bits
        if codeword(1) > 0
            bits(1) = 1;
        else
            bits(1) = 0;
        end
    elseif num_bits == 2
        %scaling to meet average power constraints and decoding the symbols
        %to bits
        if codeword(1) > 0
            bits(1) = 1;
        else
            bits(1) = 0;
        end
        if codeword(2) > 0
            bits(2) = 1;
        else
            bits(2) = 0;
        end
    elseif num_bits == 4
        %scaling to meet average power constraints and decoding the symbols
        %to bits
        power_scaled = sqrt(P/(10*rep_number));
        if codeword(1) < -2*power_scaled
            bits(1) = 0;
            bits(2) = 0;
        elseif codeword(1) < 0
            bits(1) = 0;
            bits(2) = 1;
        elseif codeword(1) < 2*power_scaled
            bits(1) = 1;
            bits(2) = 1;
        else
            bits(1) = 1;
            bits(2) = 0;
        end
        if codeword(2) < -2*power_scaled
            bits(3) = 0;
            bits(4) = 0;
        elseif codeword(2) < 0
            bits(3) = 0;
            bits(4) = 1;
        elseif codeword(2) < 2*power_scaled
            bits(3) = 1;
            bits(4) = 1;
        else
            bits(3) = 1;
            bits(4) = 0;
        end
    elseif num_bits == 6
        %scaling to meet average power constraints and decoding the symbols
        %to bits
        power_scaled = sqrt(P/(42*rep_number));
        if codeword(1) < -6*power_scaled
            bits(1) = 0;
            bits(2) = 0;
            bits(3) = 0;
        elseif codeword(1) < -4*power_scaled
            bits(1) = 0;
            bits(2) = 0;
            bits(3) = 1;
        elseif codeword(1) < -2*power_scaled
            bits(1) = 0;
            bits(2) = 1;
            bits(3) = 1;
        elseif codeword(1) < 0
            bits(1) = 0;
            bits(2) = 1;
            bits(3) = 0;
        elseif codeword(1) < 2*power_scaled
            bits(1) = 1;
            bits(2) = 1;
            bits(3) = 0;
        elseif codeword(1) < 4*power_scaled
            bits(1) = 1;
            bits(2) = 1;
            bits(3) = 1;
        elseif codeword(1) < 6*power_scaled
            bits(1) = 1;
            bits(2) = 0;
            bits(3) = 1;
        else
            bits(1) = 1;
            bits(2) = 0;
            bits(3) = 0;
        end
        
        if codeword(2) < -6*power_scaled
            bits(4) = 0;
            bits(5) = 0;
            bits(6) = 0;
        elseif codeword(2) < -4*power_scaled
            bits(4) = 0;
            bits(5) = 0;
            bits(6) = 1;
        elseif codeword(2) < -2*power_scaled
            bits(4) = 0;
            bits(5) = 1;
            bits(6) = 1;
        elseif codeword(2) < 0
            bits(4) = 0;
            bits(5) = 1;
            bits(6) = 0;
        elseif codeword(2) < 2*power_scaled
            bits(4) = 1;
            bits(5) = 1;
            bits(6) = 0;
        elseif codeword(2) < 4*power_scaled
            bits(4) = 1;
            bits(5) = 1;
            bits(6) = 1;
        elseif codeword(2) < 6*power_scaled
            bits(4) = 1;
            bits(5) = 0;
            bits(6) = 1;
        else
            bits(4) = 1;
            bits(5) = 0;
            bits(6) = 0;
        end
    end
end