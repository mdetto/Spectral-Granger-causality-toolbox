%reproduce numerical example from Wilson 1972

clear
n=200;
theta=(0:pi/(n/2):pi);
g0=[10.37 4.43; 4.43 7.32];
g1=[-3.82 2.43; -0.63 4.08];
g2=[2.07 3.24; 2.56 3.06];
g3=[0.60 0.70; -0.50 0.70];

S(1,1,:) = g0(1,1)*exp(i*0*theta)+g1(1,1)*exp(-i*1*theta)+g2(1,1)*exp(-i*2*theta)+g3(1,1)*exp(-i*3*theta)+...
                            +g1(1,1)*exp(i*1*theta)+g2(1,1)*exp(i*2*theta)+g3(1,1)*exp(i*3*theta);
S(2,1,:) = g0(2,1)*exp(i*0*theta)+g1(2,1)*exp(-i*1*theta)+g2(2,1)*exp(-i*2*theta)+g3(2,1)*exp(-i*3*theta)+...
                            +g1(1,2)*exp(i*1*theta)+g2(1,2)*exp(i*2*theta)+g3(1,2)*exp(i*3*theta);
     
S(1,2,:) = g0(1,2)*exp(i*0*theta)+g1(1,2)*exp(-i*1*theta)+g2(1,2)*exp(-i*2*theta)+g3(1,2)*exp(-i*3*theta)+...
                            +g1(2,1)*exp(i*1*theta)+g2(2,1)*exp(i*2*theta)+g3(2,1)*exp(i*3*theta);
S(2,2,:) = g0(2,2)*exp(i*0*theta)+g1(2,2)*exp(-i*1*theta)+g2(2,2)*exp(-i*2*theta)+g3(2,2)*exp(-i*3*theta)+...
                            +g1(2,2)*exp(i*1*theta)+g2(2,2)*exp(i*2*theta)+g3(2,2)*exp(i*3*theta);
                       
                                                
[SIGMA,H,B1,B2,B3]=wilson(S,20,1e-4);
     

figure(1)
csw(:,1)=H(1,2,:);
subplot(121)
plot(theta,real(csw))
set(gca,'xlim',[0 pi])
ylabel('Real W_{21}')
subplot(122)
plot(theta,imag(csw))
xlabel('angular frequency')
ylabel('Imag W_{21}')
set(gca,'xlim',[0 pi])



