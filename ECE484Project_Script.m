% ECE 484 TAKE-HOME PROJECT
% NICOLAS FAJARDO AND JOHN GENTRY

clc
clear

% Define number of elements in the array
N = 10;
if ( mod(N,2) ~= 0 )
    disp('Error. N must be even')
    return;
end

% Define the array spacing
d = 0.5;                    % NOTE: in terms of lambda (ie: 0.5 = lambda/2)

% Define PSI for the array factor calculation
PSI = 2*pi*d;

AF = zeros(1,181);



for theta = 0:1:180;
    for n = 1:N
        AF(1,theta+1) = AF(1,theta+1) + exp(1i*(n-1)*PSI*cosd(theta));
    end
end

theta = 0:1:180;

[AF_max,AF_index] = max(abs(AF));

AF_norm = AF(1,:)./AF_max(1);

figure
plot(0:1:180,abs(AF));
title('AF')

hold on

figure
plot(0:1:180,abs(AF_norm));
title('Normalized')

AF_norm_dB = 10.*log10(abs(AF_norm));

figure
plot(theta,AF_norm_dB);
title('dB')
axis([0 180 -60 0])

hold off


