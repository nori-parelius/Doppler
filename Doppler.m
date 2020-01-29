% Identifying the speed of a car using the Doppler shift in the recorded
% sound it makes when passing by. Then the approximate closest distance
% between the car and the microphone is estimated

% Loading the file
w=audioread('fysikk.wav'); 
% Preparing To Fourier transform the data
Ts=60/length(w); %sampling interval in seconds
fs = 1/Ts; %sampling frequency
% The car is heard approaching and departing twice, between 7 - 14 s and 19
% - 27 s.
w1=w(7/Ts:14/Ts); % using only data when the car is heard
w2=w(19/Ts:27/Ts);
nl=100; %number of sampling intervals
k1=floor(length(w1)/100); %number of samples in each interval
k2=floor(length(w2)/100);
N1=length(w1);
N2=length(w2);
df1=fs/N1;
df2=fs/N2;
F1=-fs/2:df1:fs/2-df1;
F2=-fs/2:df2:fs/2-df2;
t1=linspace(7,14,100);
t2=linspace(19,27,100);

%%
% Calculating the Fourier transforms of each sample
WW1=zeros(nl, N1);
WW2=zeros(nl, N2);
for i=1:nl
    WW1(i,:) = abs(fftshift(fft(w1((i-1)*k1+(1:k1)),N1)));
    WW2(i,:) = abs(fftshift(fft(w2((i-1)*k2+(1:k2)),N2)));
end

%%
% Plotting
figure(1);
subplot(2,1,1);
imagesc(F1,t1,WW1,[0,1])
colormap('gray');
title('First'), xlabel('Frequency [Hz]'), ylabel('Time [s]')
subplot(2,1,2);
imagesc(F2,t2,WW2,[0,1]);
colormap('gray');
title('Second'), xlabel('Frequency [Hz]'), ylabel('Time [s]')

%%
% Calculating the speed of the car according to the Doppler formula
% v=\frac{f_{app}-f_{dep}}{f_{app}+f_{dep}}c

fd1 = 2357;
fa1 = 2697;
fd2 = 1143;
fa2 = 1381;
c = 343*3.6; % speed of sound in km/h
v1 = (fa1-fd1)/(fa1+fd1)*c;
v2 = (fa2-fd2)/(fa2+fd2)*c;
disp(v1);
disp(v2);

%%
% Estimating the shortest distance of the observer from the car
f0=2500;
t=-2:0.1:2; %in s
cb=c/3.6; %in m/s
v1b=v1/3.6; %in m/s
v2b=v2/3.6;
figure(5);
subplot(2,1,2);
hold on;
for i=0:2
    d=i*10; %in m
    f=cb.*f0./(cb+v1b.^2.*t./(sqrt(d^2+(v1b.*t).^2)));
    plot(f, t);
    set(gca,'YDir','reverse');
end
title('First - Frequency shift for different road-observer distances (0, 10, 20 m)');
xlabel('Frequency [Hz]');ylabel('Time [s]');
grid on;
hold off;

subplot(2,1,1);
hold on;
for i=0:2
    d=i*10; %in m
    f=cb.*f0./(cb+v2b.^2.*t./(sqrt(d^2+(v2b.*t).^2)));
    plot(f, t);
    set(gca,'YDir','reverse');
end
hold off;
title('Second - Frequency shift for different road-observer distances (0, 10, 20 m)');
xlabel('Frequency [Hz]');ylabel('Time [s]');
grid on;

%%
% In the first passing of the car, the frequency shift happened during
% approximately one second, which would correspond to a distance close to
% 10 m. In the second one, the change happened over 0.8 s, which indicates
% the same distance from the observer.