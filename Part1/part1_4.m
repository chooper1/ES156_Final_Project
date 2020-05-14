%ES 156 Final Project
%Code for Part 1-4
%By: Coleman Hooper

%Parameters
%channel power constraint : P
%gaussian error standard deviation: sigma
%number of bits per transmission: num_bits
%number of simulation iterations: resolution

%E_b == P/num_bits for part 1-4

clear all;

num_bits = 6;
resolution = 1000;

%testing different P and sigma values
P_lim = [0.1;0.1;10];
sig_lim = [0.1;0.1;10];
for P=P_lim(1):P_lim(2):P_lim(3) %varying P
    for sigma = sig_lim(1):sig_lim(2):sig_lim(3) %varying sigma
        %counters to keep track of total transmissions and errors
        total = 0;
        wrong = 0;
        for j = 1:resolution
            %generate random message and encode it to symbols
            for i=1:num_bits
                bits(i) = randi([0 1]);
            end
            encoding = encode(bits, P, num_bits);
            
            %channel noise 
            n(1) = normrnd(0,sigma);
            n(2) = normrnd(0,sigma);
            h(1) = normrnd(0,1); %check this
            h(2) = normrnd(0,1); %check this
            encoding(1) = h(1)*encoding(1) +  n(1);
            encoding(2) = h(2)*encoding(2) +  n(2);
            
            %decoding the symbols to bits
            bits_out = decode(encoding, P, num_bits, h);
            %checking for bit errors
            total = total + num_bits;
            for i=1:num_bits
                if bits(i) ~= bits_out(i)
                    wrong = wrong + 1; 
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
function codeword = encode(bits, P, num_bits)
    if num_bits == 4
        %scaling to meet average power constraints and encode the message as symbols
        power_scaled = sqrt(P/10);
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
        power_scaled = sqrt(P/42);
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
function bits = decode(codeword, P, num_bits, h)
    codeword(1) = codeword(1) / h(1);
    codeword(2) = codeword(2) / h(2);
    if num_bits == 4
        %scaling to meet average power constraints and decoding the symbols
        %to bits
        power_scaled = sqrt(P/10);
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
        power_scaled = sqrt(P/42); 
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
