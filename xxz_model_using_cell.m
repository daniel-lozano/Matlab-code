sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
xx=kron(sx,sx); 
yy=kron(sy,sy); 
zz=kron(sz,sz);

N=10;% N debe ser mayor a 3

Dx=1.0;Dy=1.0;Dz=2;% anisotropia en las direcciones dadas

hx=0.0;hy=0.0;hz=0.4; %posibles direcciones y magnitudes de campo magnetico aplicado

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


%Caso de condiciones periodicas, este se adiciona cuando PBC es 1
Case_pbc=Dx*kron(sx,kron(eye(2.^(N-2)),sx)) ...;
    +Dy*kron(sy,kron(eye(2.^(N-2)),sy)) ...;
    +Dz*kron(sz,kron(eye(2.^(N-2)),sz));

Ham=Ham+PBC*(Case_pbc);%se agrega el caso de condiciones periodicas

%Caso del Hamiltoniano general
for i=(0:N-2)
    
    Id1=A{i+1};
    Id2=A{N-1-i};  
    
    
    Ham=Ham+Dx*kron(Id1,kron(xx,Id2)) ... 
        +Dy*kron(Id1,kron(yy,Id2)) ... 
        +Dz*kron(Id1,kron(zz,Id2))   ;
    
end


%Caso de campo homogeneo
if (HOMOG==1 && STAG==0 && RANDOM==0)
    disp("Caso homogeneo");
for i=(0:N-1)
    Id1=A{i+1}; 
    Id3=A{N-i};
    Ham=Ham+hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3)) ;
end


%Caso de campo alternado
elseif (HOMOG==0 && STAG==1 && RANDOM==0)
    disp("Caso Staggered");
for i=(0:N-1)
    Id1=A{i+1};  
    Id3=A{N-i};
    Ham=Ham+ ((-1.0).^(i)) *(hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3))) ;
end

%Caso de campo aleatorio
elseif (HOMOG==0 && STAG==0 && RANDOM==1)
    disp("Caso aleatorio");
    
    rng('shuffle');%cambia la semilla a una dependiente del tiempo
    
    
for i=(0:N-1)
    Id1=A{i+1};
    Id3=A{N-i};
    Ham=Ham+ (2*rand()-1) *(hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3))) ;
end

else
    disp("No es claro que caso de campo es necesario")
end


%Hallando valores propios, magnetizacion y otros valores......................



disp("Valores propios");

Vp=eig(Ham); %Halla los valores propios Vp 
[V,D]=eig(Ham); %Halla los vectores propios V y los da en una matriz
% y tambien los valores propios en una matriz diagonal D

dots=(1:dim);
plot(dots,Vp,'o')
disp(Vp(1));

GS=V(:,1);
%disp(GS);

M=zeros(dim,dim);
Corr=zeros(N,N);

for i=(0:N-1)
    Id1=A{i+1};
    Id2=A{N-i};
    M=M+kron(Id1,kron(sz,Id2));
end

Mag=GS'*M*GS;

disp("Magnetizacion="+num2str(Mag))


%Hallando la matriz de correlacion para el estado base.


for i=(0:N-1)
    Id1=A{i+1};
    Id2=A{N-i};
    
    for j=(0:N-1)
        
       Id3=A{j+1};
       Id4=A{N-j};
       
       s1=kron(Id1,kron(sz,Id2));
       s2=kron(Id3,kron(sz,Id4));
       
       Corr(i+1,j+1)= GS'*s1*s2*GS; 
       %Corr(j+1,i+1)=-GS'*s1*s2*GS;
        
    end
end
disp(Corr)
surf(Corr)





        
