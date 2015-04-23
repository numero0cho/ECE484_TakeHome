% ECE 484 TAKE-HOME PROJECT
% NICOLAS FAJARDO AND JOHN GENTRY

clc
clear

% Define number of elements in the array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 12;
if ( mod(N,2) ~= 0 )
    disp('Error. N must be even')
    return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define the array spacing
% NOTE: in terms of lambda (ie: 0.5 = lambda/2)
d = 0.5;                    


% Define PSI for the array factor calculation
PSI = 2*pi*d;


% Initialize vector for array factor data
AF = zeros(1,1801);

% Defining length of vectors
theta = 0:.1:180;
theta_deg = 0:.1:180;




% For Tschebychev Pattern, uncomment this section and comment ULA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_n = [ 1, 1.08244, 1.50839, 1.90128, 2.2031, 2.36733, 2.36733, ...
%     2.2031, 1.90128, 1.50839, 1.08244, 1 ];
% 
% for i = 1:length(theta)-1;
%     for n = 1:N
%         AF(1,i) = AF(1,i) + a_n(n)*exp(1i*(n-1)*PSI*cosd(theta_deg(i)));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% For regular ULA pattern, uncomment this section and comment Tschebychev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(theta)-1;
    for n = 1:N
        AF(1,i) = AF(1,i) + exp(1i*(n-1)*PSI*cosd(theta_deg(i)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% For random failure in 12-element Dolph-Chebyshev array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % For a completely random element failure, use the shuffled rng
% rng('shuffle')
% 
% % To test failures for all elements, manually seed the rng with a
% % different number in each instance
% % rng(11)
% 
% a = 1;
% b = 12.99;
% rand_num = floor( (b-a).*rand(4,1) + a );
% 
% % To check if the random numbers are being generated in the correct range
% % rand_num = floor( (b-a).*rand(1000,1) + a );
% % r_range = [min(rand_num) max(rand_num)];
% 
% a_n = [ 1, 1.08244, 1.50839, 1.90128, 2.2031, 2.36733, 2.36733, ...
%     2.2031, 1.90128, 1.50839, 1.08244, 1 ];
% 
% for i = 1:length(theta)-1;
%     for n = 1:N
%         if (n == rand_num(1) || n == rand_num(2) || n == rand_num(3) || n == rand_num(4))
%             AF(1,i) = AF(1,i) + 0*exp(1i*(n-1)*PSI*cosd(theta_deg(i))); 
%         else
%             AF(1,i) = AF(1,i) + a_n(n)*exp(1i*(n-1)*PSI*cosd(theta_deg(i)));
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Calculation of maximum value for array factor and its value, as well as
% calculation of normalized array factor array and normalized dB A.F.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AF_max,AF_index] = max(abs(AF));

AF_norm = AF(1,:)./AF_max;
AF_norm_dB = 20.*log10(abs(AF_norm));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% The following prints the normalized plot of the array factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot(0:1:180,abs(AF_norm));
% title('Normalized')
% hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% The following prints the dB plot of the normalized array factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(theta,AF_norm_dB);
title('dB vs \theta')
xlabel('\theta in degrees')
ylabel('dB')
axis([0 180 -40 0])
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% DEFINE NULL AND LOBE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change these values based on the graphical analysis. By changing these
% values properly, the program will be able to calculate the peaks using
% the following algorithm

nulls = N/2;                              % number of nulls on one side of peak
lobes = nulls;                            % number of lobes on one side of peak

% END NULL AND LOBE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% INTIIALIZE VALUES TO BE USED IN LOOP CALCULATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lobeData = zeros(lobes,2);              % create matrix to store lobe data
numZeros = 0;                           % initialize number of zeros measured
maxU = -100;                            % initialize maximum U measured

zeroIndex = 1;                          % initialize the index of first zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% The following loop calculates the maximum value between two nulls (or
% the start and the first null), which corresponds to a lobe. The algorithm
% finds a null my measuring if a point is the lower than the two
% neighboring points. Once a null is found, the lobeData matrix is updated
% with the measured values for the previous lobe.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:lobes
    i = zeroIndex;
    while (numZeros < nulls)
        currU = AF_norm_dB(i);
        if (currU > maxU)
            maxU = currU;
            maxIndex = i;
        end
        
        if (i > 2 && AF_norm_dB(i-1) > currU && AF_norm_dB(i+1) > currU)
           numZeros = numZeros + 1; 
           zeroIndex = i + 1;
           break;
        end 
        
        if (i == AF_index)   % if we have reached the location of absolute max
            break;
        end
        
        i = i + 1;
    end
    
    lobeData(j,1) = theta(maxIndex);
    lobeData(j,2) = maxU;
    maxU = -100;                        % reset maxU for next iteration
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Directivity calculations and outputting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_norm = (AF_norm);
U_sin = data_norm.*(sin(theta_deg.*pi./1800));
theta = 0:(pi/1800):pi;
U_s = U_sin(1:1:end);
Q = 2*pi*trapz(theta, U_s);

directivity_dim_plot1 = abs(4 * pi .* (1/Q));
directivity_db_plot1 = abs(10.*log10(directivity_dim_plot1));

fprintf('Directivity Dimensionless = %g\n', directivity_dim_plot1);
fprintf('Directivity dB = %g\n', directivity_db_plot1);
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Compute the 3 dB beamwidth from the major lobe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
halfPower_index = 0;
point = 900;
while(halfPower_index == 0)
    if(AF_norm_dB(point-1) <= -3)
       halfPower_index = 1;
       fprintf('3 dB beamwidth is %g degrees\n', (90-theta_deg(point-1))*2);
    end
    point = point - 1;
end
fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Printing of sidelobes and their peak relative to the major lobe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Angles of the sidelobes on plot 1 (90 degrees is the major lobe):');
disp('        degrees        gain in dB');
for j = 1:1:N/2
   fprintf('\t\t  %g          %g', lobeData(j,1), lobeData(j,2));
   fprintf('\n'); 
end

for j = N/2-1:-1:1
   fprintf('\t\t  %g          %g', 180-lobeData(j,1), lobeData(j,2));
   fprintf('\n'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










