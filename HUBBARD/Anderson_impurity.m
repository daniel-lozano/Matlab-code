clear; clc;

J=1;
U=2;
mu=0;
dim_site=4;%Dimension del espacio (fija)
N=5;%Numero de sitios
Np=N-1;%numero de sitios diferentes a la impureza en la red

J_v=ones(1,Np);
dim=(dim_site).^N;%dimension del espacio de Hilbert
disp("Parametros")
disp("J="+num2str(J)+", U="+num2str(U)+", mu="+num2str(mu)+", Sitios="+num2str(N)+", sitios libres="+num2str(Np))

%Llenando hopping por sitios
J_VALORES=[0.1,4.9,0.7,0.94];
for i=(1:Np)
   J_v(i)=J_VALORES(i); 
end

%Creando matrices para estados Fermionicos
C_d_do=[0,0,0,0;
        0,0,0,0;
        1,0,0,0;
        0,-1,0,0];
N_do=[0,0,0,0;
      0,0,0,0;
      0,0,1,0;
      0,0,0,1];    
    
C_d_up=[0,0,0,0;
        1,0,0,0;
        0,0,0,0;
        0,0,1,0];
N_up=[0,0,0,0;
      0,1,0,0;
      0,0,0,0;
      0,0,0,1];    
    
C_do=C_d_do';    
C_up=C_d_up'; 

P=[1,0,0,0;
   0,-1,0,0;
   0,0,-1,0;
   0,0,0,1]; 
%Creando arreglo de matrices unidad
A_id=cell(1,N);

for i=(0:N-1)
    A_id{i+1}=eye(dim_site.^i);
end

%Creando arreglo de producto kron entre matrices paridad
P_m=cell(1,Np);
P_m{1}=1;

for i=(2:Np)
P_m{i}=kron(P_m{i-1},P);    

end


disp("Comenzando el calculo")
[Ham,Num]=FERMI_HUBBARD_ANDERSON_FUNCTION(P,A_id,P_m,J_v,U,mu,dim,N,Np);
[V,D]=eig(Ham);
GS=V(:,1);
E0=D(1,1);

N_prom=GS'*Num*GS;

disp("E_0 (estado base)="+num2str(E0))
disp("N_prom (estado base)="+num2str(N_prom))

% disp("Mirando el estado base")
% b=sparse(GS);
% disp(b)
Numero_do=zeros(N,1);
Numero_up=zeros(N,1);
Numero_up_do=zeros(N,1);


for i=(1:N)
   ID1=A_id{i};
   ID2=A_id{N-i+1};
   Numero_do(i)=GS'*kron(ID1,kron(N_do,ID2))*GS; 
   Numero_up(i)=GS'*kron(ID1,kron(N_up,ID2))*GS;
   Numero_up_do(i)=GS'*kron(ID1,kron(N_up*N_do,ID2))*GS;
end
disp("Numeros de ocupacion")
disp("down")
disp(Numero_do) 
disp("up")
disp(Numero_up)
disp("Numero de doble ocupacion")
disp(Numero_up_do)





