%ES 156 Final Project
%Code for Part 3
%By: Coleman Hooper

%Parameters
%channel power constraint : P
%gaussian error standard deviation: sigma
%number of bits per transmission: num_bits
%scaling factor: rep_number

clear all;
rgb_img = imread('input.jpg');
I2 = imresize(rgb_img,[400 600]); % resize to 400x600
imwrite(I2, 'part7_initial.png');

%https://www.mathworks.com/help/matlab/creating_plots/working-with-8-bit-and-16-bit-images.html#f2-18707
I = .2989*rgb_img(:,:,1)+.5870*rgb_img(:,:,2)+.1140*rgb_img(:,:,3); %monochrome luminescence 

%800 x 1200 image, so 7680000 bits

I2 = im2double(I);
I2 = imresize(I2,[400 600]);
%apply dct
T = dctmtx(8);
B = blkproc(I2,[8 8],'P1*x*P2',T,T');
%binary mask
mask = [1   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0
          0   0   0   0   0   0   0   0];
B2 = blkproc(B,[8 8],'P1.*x',mask);

%method being used for transmission
num_bits = 4;
rep_number = 3; 
P = 10;
sigma = 1;

B2 = uint8(B2);
for k=1:size(B2,1)/8
    for l=1:size(B2,2)/8  
        %3 lowest frequency components
        c = B2((k-1)*8+1, (l-1)*8+1);
        d = B2((k-1)*8+2, (l-1)*8+1);
        e = B2((k-1)*8+1, (l-1)*8+2);
        bits = de2bi(c);
        i = length(bits)+1;
        while i <= 4
            bits(i) = 0;
            i = i + 1;
        end
        bits_temp = de2bi(d);
        i = length(bits_temp)+1;
        while i <= 4
            bits_temp(i) = 0;
            i = i + 1;
        end
        for i=1:4
            bits(i+4) = bits_temp(i); 
        end
        bits_temp = de2bi(e);
        i = length(bits_temp)+1;
        while i <= 4
            bits_temp(i) = 0;
            i = i + 1;
        end
        for i=1:4
            bits(i+8) = bits_temp(i); 
        end
        
        for i=1:12 %to get around type errors
            if bits(i) == 1
                bitvec(i) = 1;
            else
                bitvec(i) = 0;
            end 
        end
        
        %transmit the coefficients
        for i = 1:3 %repetition-3 code
            for j=1:3 %4 dct bit values
                %encoding the message to symbols
                encoding = encode_channel(bitvec((i-1)*4+1:i*4), P, num_bits, rep_number);
                %channel noise
                n(1) = normrnd(0,sigma);
                n(2) = normrnd(0,sigma);
                h(1) = normrnd(0,1); 
                h(2) = normrnd(0,1); 
                encoding(1) = h(1)*encoding(1) +  n(1);
                encoding(2) = h(2)*encoding(2) +  n(2);
                %decoding the symbols to bits
                bits_out_temp(i,(j-1)*4+1:j*4) = decode_channel(encoding, P, num_bits, rep_number, h);
            end
        end
        %recovering the repeated transmissions and using majority rule
        for j=1:3
            bits_out = bits_out_temp(j,:);
            for i=1:4
                count = 0;
                if(bits_out(i) == 1)
                    count = count + 1;
                end
                if(bits_out(i+4) == 1)
                    count = count + 1;
                end
                if(bits_out(i+8) == 1)
                    count = count + 1;
                end
                if count > 1
                    bits_maj((j-1)*4+i) = 1;
                else
                    bits_maj((j-1)*4+i) = 0;
                end
            end
        end
        %set 8 by 8 block using the recovered frequency components
        B3((k-1)*8+1:8*k,(l-1)*8+1:8*l) = zeros(8);
        int_out = bi2de(bits_maj(1:4));
        B3((k-1)*8+1,(l-1)*8+1) = int_out;
        int_out = bi2de(bits_maj(5:8));
        B3((k-1)*8+2,(l-1)*8+1) = int_out;
        int_out = bi2de(bits_maj(9:12));
        B3((k-1)*8+1,(l-1)*8+2) = int_out;
    end
end

%apply the inverse DCT
I3 = blkproc(B3,[8 8],'P1*x*P2',T',T);
%show the original image and the recovered image
imshow(I2), figure, imshow(I3);

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
