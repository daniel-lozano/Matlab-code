sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
xx=kron(sx,sx); 
yy=kron(sy,sy); 
zz=kron(sz,sz);
N=8;% N debe ser mayor a 3

Dx=1.0;Dy=1.0;Dz=2.0;% anisotropia en las direcciones dadas

hx=0.0;hy=0.0;hz=0.0; %posibles direcciones y magnitudes de campo magnetico aplicado

dim=2.^N; %dimension del espacio a analizar

%AGREGANDO PARAMETROS DE CAMPO MAGNETICO, CONDICIONES DE FRONTERA
HOMOG=0 ;
STAG=1;
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
t=0;
Ham=Ham+Hamiltonian_t(N,dim,Dx,Dy,Dz,hx,hy,hz,A,t,HOMOG,STAG,RANDOM,PBC);

%Hallando valores propios, magnetizacion y otros valores......................

Vp=eig(Ham); %Halla los valores propios Vp 
[V,D]=eig(Ham); %Halla los vectores propios V y los da en una matriz
% y tambien los valores propios en una matriz diagonal D

GS=V(:,1);%Vector de estado base

%evolucionando!
disp("Realizando evolucion temporal........................................");
disp("Parametros");
N_T=200;
disp("Numero de pasos="+num2str(N_T));

T_i=0;
disp("Tiempo inicial="+num2str(T_i));

T_f=2;
disp("Tiempo final="+num2str(T_f));

Time=linspace(T_i,T_f,N_T);
dT=Time(2)-Time(1);
disp("Paso="+num2str(dT));

GS_t=GS;

disp("Estado base en t=0");
disp(GS_t);
disp("Norma="+num2str(GS_t'*GS_t));

for i=(1:N_T-1)
    old=GS_t;%copia el ultimo valor de la evolucion
    
    t=Time(i)+dT/2.0;%tiempo en el que sera evolucionado
    
    Matrix=Hamiltonian_t(N,dim,Dx,Dy,Dz,hx,hy,hz,A,t,HOMOG,STAG,RANDOM,PBC);%Calcula matriz de evolucion en tiempo dado
    
    EXP=expm(-1i*Matrix*dT);% e^{-iH(t_i+Delta_t/2)\Delta t}
    
    GS_t=EXP*old;%evolucionando
    
    GS_t=GS_t/sqrt(GS_t'*GS_t);%normalizando
end

%Normalizando el resultado
disp("Estado base en t=T_f")
disp(GS_t)
disp("Norma="+num2str(GS_t'*GS_t));

%Hallando la fidelidad
disp("Hallando la fidelidad del estado final con el estado inicial")
fid=GS_t'*GS;
disp("Valor complejo")
disp(fid)
disp("Normal")
disp(fid'*fid)

%Estado aleatorio..........................................................
up=[1,0];
down=[0,1];

a_new=up;
%Construyendo estado aleatorio
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

Total_time=2;
N_T=200;% da 0.45
Time=linspace(0,Total_time,N_T);
dT=Time(2)-Time(1);
r_s_t=random_state;

%evolucionando el estado aleatorio escogido

for i=(1:N_T-1)
    old=r_s_t;%copia el ultimo valor de la evolucion
    
    t=0%Time(i)+dT/2.0;%tiempo en el que sera evolucionado
    
    Matrix=Hamiltonian_t(N,dim,Dx,Dy,Dz,hx,hy,hz,A,t,HOMOG,STAG,RANDOM,PBC);%Calcula matriz de evolucion en tiempo dado
    
    EXP=expm(Matrix*dT);% e^{-iH(t_i+Delta_t/2)\Delta t}
    
    r_s_t=EXP*old;%evolucionando
    
    r_s_t=r_s_t/sqrt(r_s_t'*r_s_t);%normalizando
end

%Hallando la fidelidad
disp("Hallando la fidelidad del estado final aleatorio con el estado inicial")
fid=r_s_t'*GS;
disp("Normal")
disp(fid'*fid)









        
