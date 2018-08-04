sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
xx=kron(sx,sx); 
yy=kron(sy,sy); 
zz=kron(sz,sz);
N=4;% N debe ser mayor a 3

Dx=1.0;Dy=1.0;Dz=-5.0;% anisotropia en las direcciones dadas

hx=0.0;hy=0.0;hz=-1.0; %posibles direcciones y magnitudes de campo magnetico aplicado

dim=2.^N; %dimension del espacio a analizar

%AGREGANDO PARAMETROS DE CAMPO MAGNETICO, CONDICIONES DE FRONTERA
HOMOG=1 ;
STAG=0;
RANDOM=0;
Realizations=10;

PBC=0; %Periodic Boundary Conditions

Ham=zeros(dim,dim);

%Llenando una celda con las posibles matrices identidad necesarias 
A=cell(1,N);
for i=(0:N-1)
    A{i+1}=eye(2.^i);
end


%Funcion que da el Hamiltoniano

Ham=Ham+Hamiltonian(N,dim,Dx,Dy,Dz,hx,hy,hz,A,HOMOG,STAG,RANDOM,PBC);

%Hallando valores propios, magnetizacion y otros valores......................

Vp=eig(Ham); %Halla los valores propios Vp 
[V,D]=eig(Ham); %Halla los vectores propios V y los da en una matriz
% y tambien los valores propios en una matriz diagonal D

GS=V(:,1);%Vector de estado base

%Estado aleatorio..........................................................
up=[1,0];
down=[0,1];

a_new=up;
%Construyendo estado antiferromagnetico
for i=(2:N)
    
    if(mod(i,2)==0)
        step=down;
    
    else
        step=up;
    end
    
    
    a_new=kron(a_new,step);
end
random_state=a_new';
%disp(size(random_state))
%disp(random_state)

Total_time=15;
N_T=3000;%
Time=linspace(0,Total_time,N_T);
dT=Time(2)-Time(1);
r_s_t=random_state;
    
EXP=expm(Ham*dT); % e^{-iH(t_i+Delta_t/2)\Delta t

%evolucionando el estado aleatorio escogido

for i=(1:N_T-1)
    old=r_s_t;%copia el ultimo valor de la evolucion
    
    r_s_t=EXP*old;%evolucionando
    
end

r_s_t=r_s_t/sqrt(r_s_t'*r_s_t);%normalizando

%Hallando la fidelidad
disp("Hallando la fidelidad del estado final aleatorio con el estado inicial")
fid=r_s_t'*GS;
disp("Parte real")
disp(real(fid))
disp(r_s_t)
disp(GS)







        
