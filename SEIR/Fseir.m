

function dY= Fseir(t,Y,alpha,beta,gamma,delta,mu,N,Nh,p,Nv,Sv,Ev,V)

S=Y(1);
E=Y(2);
I=Y(3);
R=Y(4);

dY = zeros(2,1);

dY(1) = alpha*(1-p)*Nh-((beta*S*(Nv-Sv-Ev))/Nh)-mu*S;
dY(2) =  ((beta*S*(Nv-Sv-Ev))/Nh)-delta*E-mu*E;
dY(3) =  delta*E-gamma*I-mu*I;
dY(4)=   alpha*V*N+gamma*I-mu*R; % its not necessary this line (dR+dI+dE+dS=0)

end