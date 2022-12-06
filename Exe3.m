% 3.a

T=10;
t=0:1/100:T;
n1=-20:1:20;
n2=0:1:19;
phi = exp(1i.*(2*pi/T)*(t').*n1);
psy_matrix=zeros(length(t),20);
for i=0:19                                                  
    psy_matrix(:,i+1)=rectangularPulse(T*i/20,T*(i+1)/20,t);
end
ft = 4*cos(((4*pi)/T).*t) + sin(((10*pi)/T).*t);
gt = 2*sign(sin(((6*pi)/T).*t)) - 4*sign(sin(((4*pi)/T).*t));
%projections
cn_1=projection_res(ft,phi,T);
table(round(cn_1))
cn_2=projection_res(ft,psy_matrix,T);
table(cn_2)
Cn_3=projection_res(gt,phi,T);
table(Cn_3)
Cn_4=projection_res(gt,psy_matrix,T);
table(Cn_4)


%build f on psy: 
f_hat1=zeros(length(t),1);
for i=1:41
    f_hat1=f_hat1+cn_1(i).*phi(:,i);     
end
f_hat2=zeros(length(t),1);
for i=1:20
    f_hat2=f_hat2+cn_2(i).*psy_matrix(:,i);     
end

plot(t,f_hat2,'k',Linewidth=2);
hold on
xlabel('t [sec]')
ylabel('Amp[v]')
title('f(t) based on \psi_n(t) ')
plot(t,ft,'--',LineWidth=3);
legend([{'f(t)'};{'f_build'}]);
grid minor
hold off 
%build f on phi 
plot(t,f_hat1,'k',Linewidth=2);
hold on
xlabel('t [sec]')
ylabel('Amp[v]')
title('f(t) based on \phi_n(t) ')
plot(t,ft,'--',LineWidth=3);
legend([{'f_build'};{'f(t)'}]);
grid minor
hold off 

%build g on psy: 
g_hat1=zeros(length(t),1);
for i=1:41
    g_hat1=g_hat1+Cn_3(i).*phi(:,i);     
end
g_hat2=zeros(length(t),1);
for i=1:20
    g_hat2=g_hat2+Cn_4(i).*psy_matrix(:,i);     
end
plot(t,g_hat2,'--',Linewidth=2);
hold on
xlabel('t [sec]')
ylabel('Amp[v]')
title('g(t) based on \psi_n(t) ')
plot(t,gt,'k',LineWidth=1);
legend([{'g_build'};{'g(t)'}]);
grid minor
hold off 

%build f on phi 
plot(t,g_hat1,'--',Linewidth=2);
hold on
xlabel('t [sec]')
ylabel('Amp[v]')
title('g(t) based on \phi_n(t) ')
plot(t,gt,'k',LineWidth=1);
legend([{'g_build'};{'g(t)'}]);
grid minor
hold off 


function [Coeff] = projection_res(x_sample_T,matrix_phi,T)
    t1=0:0.01:T;
    Coeff=zeros(size(matrix_phi,2),1);
    for i=1:size(matrix_phi,2)
        denom=trapz(t1,(matrix_phi(:,i)).*conj(matrix_phi(:,i)));
        numer=trapz(t1,x_sample_T'.*conj(matrix_phi(:,i)));
     
        Coeff(i)=numer/denom;
    end
end


