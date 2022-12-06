% 1.a
wm=3.*pi;
ti=0.2:0.01:3;
xt=(4./(wm*pi*ti.^2)).*((sin(wm*ti)).^2).*cos(wm*ti).*sin(2*wm*ti);
plot(ti,abs(xt),LineWidth=2);
hold on
xlabel('time[sec])','FontSize',12);
ylabel('|X(t)|[V]','FontSize',12);
title('Exe.1.a - X(t)')
grid minor
hold off

% 1.b
wk=(-17*pi:0.1*pi:17*pi);
%  after trigonometry the new xt:  
%  xt_new=((4.*wm)/pi).*(sinc((wm.*ti)/pi).^2).*(1/4.*1i)*(exp(1i.*3.*wm.*ti)-exp(-1i.*3.*wm.*ti)+exp(1i.*wm.*ti)-exp(1i.*wm.*ti))
%  in stages we will abstract the fourier transform of the signal
xfw1=triangularPulse(wm,3*wm,5*wm,wk);
xfw2=-triangularPulse(-5*wm,-3*wm,-wm,wk);
%  this two signals centered in -3wm and 3wm
xfw3=+triangularPulse(-wm,wm,3*wm,wk);
xfw4=-triangularPulse(-3*wm,-wm,wm,wk);
%  this two signals centered in -wm and wm 
plot(wk,xfw1,LineWidth=1)
hold on
plot(wk,xfw2,LineWidth=1)
plot(wk,xfw3,LineWidth=1)
plot(wk,xfw4,LineWidth=1)
grid minor
title('Exe.1.b.1')
legend('3wm','-3wm','-wm','wm')
xlabel('w[rad/sec]','FontSize',12);
ylabel('|Trianglef(w)|','FontSize',12);
hold off

%  now lets do the summary of all these signals and we get: 
%  in stages we will abstract the fourier transform of the signal
plot(wk,abs(xfw1+xfw2+xfw3+xfw4),LineWidth=1)
grid minor
title('Exe 1.b.2 ')
xlabel('w[rad/sec]','FontSize',12);
ylabel('|Xf(w)|','FontSize',12);
hold off


% 1.c.1
%  the maximal frequency is for this signal according to the graph is
%  (46.9219~)/(pi)=15 
%  thus the maximal frequency[w] is 15*pi. 


%  we need to suggest Ts so we can restore the signal. 
%  ts=1/f according to do that we need to understand Nyquist–Shannon term which
%  means ws=2wmax , ws=2*15*pi=30*pi 
%  because the x fourier transform had no spaces in it and the maximal
%  frequency is 15*pi so we need to sample according to Nyquist–Shannon 
%  term for actual signal without any errors during the procces
%  so if we calculate in terms of time : ts=1/fs=2*pi/(2*pi*15)=1/15 [sec]

% 1.c.2 
%  now we will sample the signal. in terms of ZOH according to the system
%  structure.
Ts_40n=1/15;
ts_vec=0.2:Ts_40n:3;
sample_vec=(4./(wm*pi*ts_vec.^2)).*((sin(wm*ts_vec)).^2).*cos(wm*ts_vec).*sin(2*wm*ts_vec);
plot(ti,xt,LineWidth=2)
hold on
plot(ts_vec,sample_vec,'*')
stem(ts_vec,sample_vec,'b'); 
xlabel('t [sec]','FontSize',12);
ylabel('X(t)[V]','FontSize',12);
title('Exe 1.c.2 - X(t) & Xs(t)')
legend([{'X(t)'};{'X_s(t)'};{'X_s(t)- using stem'}]);
grid minor
hold off 
%1.c.3
zoh_f=sample_vec(1)*rectangularPulse(ts_vec(1),ts_vec(2),ti);
for i=2:length(sample_vec)
zoh_f=zoh_f+sample_vec(i)*rectangularPulse(ts_vec(i),ts_vec(i)+Ts_40n,ti);
end 
zoh_f;
 xt=4./(wm*pi*ti.^2).*((sin(wm.*ti)).^2).*cos(wm.*ti).*sin(2.*wm.*ti);
plot(ti,xt,LineWidth=2);
hold on
plot(ti,zoh_f,LineWidth=1)
grid minor
xlabel('t[sec]','FontSize',12);
ylabel('X(t)[V]','FontSize',12);
title('Exe 1.c.3 - Xzoh(t) & X(t)')
hold off
legend('X(t)','Xzoh(t)')
%1.d 
%we need describe the X_fzoh(w) as a fourier transform of the signal
%xzoh(t). the analitical expression is the triangular signal in the
%frequency domain multiplied by e^(-j*w*Ts/2)*sinc(w/Ts*(1/2*pi))
ws=30*pi;
triangle_total=xfw1+xfw2+xfw3+xfw4;
xwzoh=(triangle_total-triangularPulse(5*wm,7*wm,7*wm,wk)+triangularPulse(-7*wm,-7*wm,-5*wm,wk)).*exp(-1i*Ts_40n*wk/2).*sinc(wk/ws);
plot(wk,abs(xwzoh),LineWidth=1)
title('Exe.1.d - XFzoh(w)')
xlabel('w[rad/sec]','FontSize',12);
ylabel('abs(X_zoh(w))');
grid minor

%1.e 
%the spectrum of the ideal reconstructor:
Hw=rectangularPulse(-ws/2,ws/2,wk).*(((exp(pi*wk*1j/ws)))./(sinc(wk/ws)));
plot(wk,abs(Hw),LineWidth=1)
hold on 
xlabel('w[rad/sec])','FontSize',12);
ylabel('abs(H(w))','Fontsize',12);
grid minor
title('Exe.1.e.1 - The spectrum of the ideal reconstructor')
hold off

%very important to see how the filter looks like. 
%Now we need to calculate the xfrec so we need to use trapz function. how
%am i suppose to do that? 
%lets multiply the filter with xwfzoh
xWrec=Hw.*xwzoh;
%now we need to do inverse fourier transform using the trapz function 
X_rec=zeros(1,length(ti));
for i=1:length(ti)
    X_rec(i)=1/(2*pi*1i)*trapz(wk,xWrec.*exp(1i*wk*ti(i)));
end
plot(ti,xt,LineWidth=2);
hold on
plot(ti,X_rec,'--',LineWidth=2)
title('Exe 1.e.2 - X(t) & X-rec(t) ');
xlabel('t[sec]','FontSize',12);
ylabel('x(t)[V]','FontSize',12);
legend([{'X-rec(t)'};{'X(t)'};]);
grid minor 
hold off
%As we can see (only visualy) the reconstructor had completed it's task
%property. we can estimate the error of it by doing some calculations. 
% 1.f 

%we saw before that :  the maximal frequency is for this signal according to the graph of exe 1.b is
%  (46.9219~)/(pi)=15 
%  thus the maximal frequency[w] is 15*pi and we have demonstratred how we could reconstruct the signal according to nyquist term with Ts=1/15
% Now we should see that we are not working in shannon nyquist frequency: 
%Ts_new=2pi/9wm=2/27 > Ts_old ==> and ws_new=27*pi (2wm>ws)=(30*pi>27*pi)
%according to nykquist term we cant reconstruct the signal.
% we will try to demonstrate the same proccess with the new ws and we will see the results . 
Ts_new=2/27;
ts_vec2=0.2:Ts_new:3;
sample_vec_new=(4./(wm*pi*ts_vec2.^2)).*((sin(wm*ts_vec2)).^2).*cos(wm*ts_vec2).*sin(2*wm*ts_vec2);
plot(ti,xt,LineWidth=2)
hold on
plot(ts_vec2,sample_vec_new,'*',LineWidth=2)
xlabel('t[sec]','FontSize',12);
ylabel('x(t)[V]','FontSize',12);
title('Exe 1.f.1')
legend([{'x(t)'};{'x-s-new(t)'}]);
grid minor
hold off 

ws_new=27*pi;
triangle_total=xfw1+xfw2+xfw3+xfw4;
xwzoh_new=(triangle_total-triangularPulse(5*wm,7*wm,7*wm,wk)+triangularPulse(-7*wm,-7*wm,-5*wm,wk)).*exp(-1i*Ts_new*wk/2).*sinc(wk/ws_new);
plot(wk,abs(xwzoh_new));
hold on
plot(wk,abs(xwzoh));
legend('xw-zoh-new(ws=27*pi)','xw-zoh(ws=30*pi)')
xlabel('w[rad/sec]','FontSize',12);
ylabel('abs(X-zoh-new(w))');
title('Exe 1.f.2 - the diffrance between the zoh functions in the frequency domain ')
grid minor
hold off

Hw_new=rectangularPulse(-ws_new/2,ws_new/2,wk).*(((exp(pi*wk*1j/ws_new)))./(sinc(wk/ws_new)));
xWrec_new=Hw_new.*xwzoh_new;
%now we need to do inverse fourier transform using the trapz function 
X_rec_new=zeros(1,length(ti));
for i=1:length(ti)
    X_rec_new(i)=1/(2*pi*1i)*trapz(wk,xWrec_new.*exp(1i*wk*ti(i)));
end
plot(ti,xt,LineWidth=2);
hold on
plot(ti,X_rec_new,'r--',LineWidth=2)
legend('x(t)','X-rec-new(t)')
title('Exe 1.f.3 - the new reconstructed signal')
xlabel('t[sec]','FontSize',12);
ylabel('x(t)[V]','FontSize',12);
grid minor
hold off
%we could not reconstruct the signal accuratly because of nyquist term .
%2_a

tl = 0:0.01:2 ;                        %continious time vector
%we found the signal period as 2sec                             
Tk = 2/15;                                                          
Wa = 7*pi;                            
Wb = 4*pi;
tsample = 0:Tk:(2-Tk);                   
Xt = 5*cos(Wa*tl) - 3*sin(Wb*tl);
Xs = 5*cos(Wa*tsample) - 3*sin(Wb*tsample);
plot(tl,Xt,Linewidth=2);
hold on;
stem(tsample,Xs,'r');
xlabel('t[sec])');
ylabel('Xi(t)[V]');
title('Exe 2.a - Continious signal X(t) & Sampeled signal Xs(t)');
legend('X(t)','Xs(t)')
grid minor
hold off

% 2.b
%in this section we need to calculate the fourier matrix of the signal 
l = -7:7;
exp_vec=zeros(1,length(l));
for k = -7:7
%creating exponent vector [1x15]
exp_vec(k+8) = exp(1i*pi*k);                
end
%creating exponent matrix [15x15]
exp_mat = exp((1i*pi*l).*tsample');              
fourier_coefficient = inv(exp_mat)*Xs.';                       
table(fourier_coefficient);

% 2.c
%creating 201X15 exponents matrix                                                
FMat = exp((1i*pi*l).*tl.');
%reconstructing from Fourier Coefficients
x_r =  FMat*fourier_coefficient;                  
plot(tl, Xt,'k',Linewidth=4)
hold on
plot(tl,x_r,'m',LineWidth=0.8)
title('Exe 2.c - reconstructed singal from fourier coefficients & orginial signal')
legend('X(t)', 'Xrec(t)')
ylabel('Amp[V]');
xlabel('t[sec]');
grid minor
hold off


% 2.d.1
N=15;
%picking random time samples

tn= sort(2*rand([N,1]));
Xn = 5*cos(Wa*tn)-3*sin(Wb*tn);
plot(tl,Xt,'r',LineWidth=2);
hold on;
stem(tn,Xn,'b',LineWidth=1.5);
title('Exe 2.d.1  - Xt signal & randomly sampeled signal - Xn');
ylabel('Amp');
xlabel('t[sec]');
legend('signal','sampeled signal')
grid minor
hold off
%the values of x(random time vector)=samples vector 
Xn;

% 2.d.2
FF_mat = zeros(15,15);
W0 = 2*pi/2;
for n = 1:15
    for k = -7:7
        FF_mat(n,k+8) = exp(1i*tn(n)*W0*k);
    end
end

% fourier coefficients vector from random sampling
coaff_rand = (inv(FF_mat))*Xn; 
table(coaff_rand);


% 2.d.3 

    X_rec = zeros(1,length(tl));
for j = 1:length(tl)
    for m = 1:length(tn)
        X_rec(j) = X_rec(j) + coaff_rand(m)*exp(1i*W0*tl(j)*(m-8));
    end
end

plot(tl,X_rec,'b','LineWidth',4); 
hold on
plot(tl,Xt,'g',LineWidth=1.5)
title('Exe 2.d.3 - signal & reconstructed signal from random samples');
ylabel('Amp[V]');
xlabel('t[sec]');
legend('Xrec','Xt');
grid minor
hold off  ;
%now we repeating from a to d with uncertainty in the samples' location
% 2.e.1.random samples
tn1 = zeros(N,1);
for n = 1:15
    tn1(n) = tn(n) + 0.01*rand(1); %t(n) random
end
tn1 = sort(tn1);
Xn_new = 5*cos(Wa*tn1)-3*sin(Wb*tn1);
plot(tl,Xt,'r',LineWidth=1.5)
hold on;
stem(tn1,Xn_new,'b');
title('Exe 2.e.1 - Xt signal & new randomly sampeled signal Xn-new(random with uncertenty)');
ylabel('Amp[V]');
xlabel('t[sec]');
grid minor
legend('signal','sampeled');
hold off

% 2.e.1.uniform samples
%time vector with fixed intervals and random intervals added
tn2_new = zeros(N,1);
for n = 1:15
    tn2_new(n) = tsample(n) + 0.01*rand(1);
end
tn2_new = sort(tn2_new);
Xuniform_new = 5*cos(Wa*tn2_new)-3*sin(Wb*tn2_new);
plot(tl,Xt,'b',LineWidth=2.5);
hold on;
stem(tn2_new,Xuniform_new,'r',LineWidth=1.5)
title('Exe 2.e.2 - Fixed-interval sampling Xt  & new signal (uniform with uncertenty)');
ylabel('Amp[V]');
xlabel('t[sec]');
grid minor
legend('signal','sampeled');
hold off

% 2.e.3.1 uniform with uncertenty matrix caculate
FFn2_mat = zeros(15,15);
W0 = 2*pi/2;
for n = 1:15
    for k = -7:7
        FFn2_mat(n,k+8) = exp(1i*tn2_new(n)*W0*k);
    end
end
fourier_coefficients_uniform = (inv(FFn2_mat))*Xuniform_new;  %coefficients vector 
four_co_uniform_uncertenty=table(fourier_coefficients_uniform)

% 2.e.3.2 random with uncertenty matrix caculate
FFn3_mat = zeros(15,15);
W0 = 2*pi/2;
for n = 1:15
    for k = -7:7
        FFn3_mat(n,k+8) = exp(1i*tn1(n)*W0*k);
    end
end
fourier_coefficients_random = (inv(FFn3_mat'*FFn3_mat))*(FFn3_mat')*Xn_new;  %coefficients vector 
four_co_random_uncertenty=table(fourier_coefficients_random)

X_recn_uniform = zeros(1,length(tl));
for j = 1:length(tl)
    for m = 1:length(tn2_new)
        X_recn_uniform(j) = X_recn_uniform(j) + fourier_coefficients_uniform(m)*exp(1i*W0*tl(j)*(m-8));
    end
end



   
X_recn_random = zeros(1,length(tl));
for j = 1:length(tl)
    for m = 1:length(tn)
        X_recn_random(j) = X_recn_random(j) + fourier_coefficients_random(m)*exp(1i*W0*tl(j)*(m-8));
    end
end

plot(tl,X_recn_uniform,'g',LineWidth=3.5)
hold on 
plot(tl,Xt,'--',LineWidth=1.5)
title('Exe 2.e.3 - original and reconstructed signal (uniform with uncertainty)');
ylabel('Amp[V]');
xlabel('t[sec]');
grid minor
legend('reconstructed','original');
hold off


plot(tl,X_recn_random,'b',LineWidth=4.5); 
hold on
plot(tl,Xt,'--',LineWidth=2)
title('Exe 2.e.4 - original and reconstructed signal (random with uncertenty');
xlabel('t[sec]')
grid minor
ylabel('Amp[V]')
legend('Xrec','Xorginial')
hold off  


% 2.e.5
%condition number---%
uniform_sampling_condition_number = cond(FFn2_mat)
random_samping_condition_number = cond(FFn3_mat)


% 2.f 
T=2;
Ts_40n = 2/40;
tn = sort(1.999*rand(40,1));
tn_40P = zeros(40,1);
for n = 1:40
    tn_40P(n) = tn(n) + 0.01*rand(1);
end
tn_40P = sort(tn_40P);
FF40 = zeros(40,15);
W0 = 2*pi/T;
for n = 1:40
    for k = -7:7
        FF40(n,k+8) = exp(1i*tn_40P(n)*W0*k);
    end
end
X40 = 5*cos(Wa*tn_40P)-3*sin(Wb*tn_40P);
a40 = (inv(FF40'*FF40))*(FF40')*X40; 
%coefficients vector 
X_rec40 = zeros(1,length(tl));
for n = 1:length(tl)
    for m = 1:15
        X_rec40(n) = X_rec40(n) + a40(m)*exp(1i*W0*tl(n)*(m-8));
    end
end
plot(tl,X_rec40 ,'b',LineWidth=3);
hold on 
plot(tl,Xt ,'--' ,LineWidth=1.5);
grid minor
title('Exe 2.f.1 reconstructed signal from 40 random samples')
xlabel('t[sec]')
ylabel('Amp[V]')
legend('reconstructed','original')
hold off


plot(tl,Xt,'r',LineWidth=2)
hold on
stem(tn_40P,X40,'b')
title('Exe 2.f.2 singal and 40 random sample of the signal')
xlabel('t[sec]')
grid minor 
ylabel('Amp[V]')
legend('signal','sampels')
hold off
cond40 = cond(FF40)


