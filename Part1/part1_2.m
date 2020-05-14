%ES 156 Final Project
%Code for Part 1-2 
%By: Coleman Hooper

%Parameters
%channel power constraint : P
%gaussian error standard deviation: sigma

%E_b == P/2 for part 1-2

clear all;

%Defining the gaussian function G(x)
fun = @(x) exp(-x.^2);

%testing different P and sigma values
P_lim = [0.1;0.1;10];
sig_lim = [0.1;0.1;10];
for P=P_lim(1):P_lim(2):P_lim(3) %varying P
    for sigma = sig_lim(1):sig_lim(2):sig_lim(3) %varying sigma
        x = sqrt(P/2)/(sqrt(2)*sigma);
        %error function complement integral
        P_e(round(P*10),round(sigma*10)) = log10((1/sqrt(pi))*integral(fun, x, Inf));
        EPB_noise(round(P*10),round(sigma*10)) = 20*log10((P/2)/(sigma^2));
    end
end

%plotting BER versus Noise Ratio (Analytically)
%note that the simulated BER versus Noise Ratio plots for 4-QAM were done
%using part2_1.m with repeat = 1.
figure;
hold on;
xlabel("Energy-Per-Bit over Noise Ratio (dB)");
ylabel("Bit Error Probability");

for i=1:1:size(P_e,1)
    scatter(EPB_noise(i,:), P_e(i,:), 1, 'k');
end
