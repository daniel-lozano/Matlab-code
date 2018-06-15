sx=[0 1;1 0]; sy=[0 -1i;1i 0]; sz=[1 0;0 -1];
xx=kron(sx,sx); 
yy=kron(sy,sy); 
zz=kron(sz,sz);

N=4;% N debe ser mayor o igual a 3

Dx=1.0;Dy=1.0;Dz=-5.5;% anisotropia en las direcciones dadas

hx=0.0;hy=0.0;hz=10.0; %posibles direcciones y magnitudes de campo magnetico aplicado

dim=2.^N; %dimension del espacio a analizar

%AGREGANDO PARAMETROS DE CAMPO MAGNETICO, CONDICIONES DE FRONTERA
HOMOG=1 ;
STAG=0;
RANDOM=0;


PBC=0; %Periodic Boundary Conditions

Ham=zeros(dim,dim);

%Caso de condiciones periodicas, este se adiciona cuando PBC es 1
Case_pbc=Dx*kron(sx,kron(eye(2.^(N-2)),sx)) ...;
    +Dy*kron(sy,kron(eye(2.^(N-2)),sy)) ...;
    +Dz*kron(sz,kron(eye(2.^(N-2)),sz));

Ham=Ham+PBC*(Case_pbc);


if(N==3)
    disp("Solucion para 3 sitios")
end

%Caso del Hamiltoniano general
for i=(0:N-2)
    
    Id1=eye(2.^(i));
    Id2=eye(2.^(N-i-2));  
    Id3=eye(2.^(N-i-1));
    
    Ham=Ham+Dx*kron(Id1,kron(xx,Id2)) ... 
        +Dy*kron(Id1,kron(yy,Id2)) ... 
        +Dz*kron(Id1,kron(zz,Id2))   ;
    
end


%Caso de campo homogeneo
if (HOMOG==1 && STAG==0 && RANDOM==0)
    disp("Caso homogeneo");
for i=(0:N-1)
    Ham=Ham+hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3)) ;
end


%Caso de campo alternado
elseif (HOMOG==0 && STAG==1 && RANDOM==0)
    disp("Caso Staggered");
for i=(0:N-1)
    Ham=Ham+ ((-1.0).^(i)) *(hx*kron(Id1,kron(sx,Id3)) ... 
        +hy*kron(Id1,kron(sy,Id3)) ... 
        +hz*kron(Id1,kron(sz,Id3))) ;
end

%Caso de campo aleatorio
elseif (HOMOG==0 && STAG==0 && RANDOM==1)
    disp("Caso aleatorio");
for i=(0:N-1)
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
disp(Vp);
GS=V(:,1);

M=zeros(dim,dim);

for i=(0:N-1)
    Id1=eye(2.^(i));
    Id2=eye(2.^(N-1-i));
    M=M+kron(Id1,kron(sz,Id2));
end

Mag=GS'*M*GS;

disp(Mag)

    



        
