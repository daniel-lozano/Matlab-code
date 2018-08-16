J=1;
U=2;
mu=0;

N_max=5;%Numero maximo de estados por sitio
sites=4;%Numero de sitios debe ser mayor a 1
PBC=0;%condiciones periodicas
dim_site=(N_max+1);%dimension de un sitio


dim=dim_site.^sites;% dimension de todo el espacio
A_id=cell(1,sites);

for i=(0:sites-1)
    A_id{i+1}=eye(dim_site.^i);%disp("tamano de A")%disp(size(A_id{i+1}))
end


A_dag=diag(sqrt(1:N_max),-1);
A=A_dag';
N=A_dag*A;


%
% disp("Op. creacion")
% disp(A_dag)
%disp("Op. destruccion")
%disp(A)
%disp("Op. Numero")
%disp(N)
%
%En la funcion se crea la matriz que representa el hamiltoniano y de igual
%manera la funcion que representa al numero total N
[Ham,Num]=BOSE_HUBBARD_FUNCTION(A_dag,A,N,A_id,J,U,0,dim,dim_site,sites,PBC);
%mirando las dimensiones
disp("dimension del espacio "+ num2str(dim))
disp("dimension del Hamiltoniano "+ num2str(size(Ham)))

[V,D]=eig(Ham);%Hallando los vectores y valores propios

GS=V(:,1);%Estado base
Eo=D(1,1);
disp("Energia del estado base="+num2str(Eo))
N_prom=GS'*Num*GS;
disp("Numero promedio de particulas="+num2str(N_prom))

for i=(1:sites)
   ID1=A_id{i};
   ID2=A_id{sites-i+1};
   disp("sitio " + num2str(i))
   Numero=kron(ID1,kron(N,ID2));
   disp(GS'*Numero*GS)
end



MU=linspace(-mu,mu,100);%posibles valores del potencial quimico
NU_MU=zeros(size(MU));%Aqui seran guardados los valores del numero promedio

for i=(1:100)
    [Ham,Num]=BOSE_HUBBARD_FUNCTION(A_dag,A,N,A_id,J,U,MU(i),dim,dim_site,sites,PBC);

    [V,D]=eig(Ham);
    GS=V(:,1);%Estado base
    NU_MU(i)=GS'*Num*GS;
end
figure(1)
plot(MU,NU_MU)
ylabel("Numero")
xlabel("mu")
hold on





