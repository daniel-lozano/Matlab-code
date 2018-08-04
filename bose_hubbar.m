J=1;
U=1;
mu=1;

N_max=1;%Numero maximo de estados por sitio
sites=4;%Numero de sitios debe ser mayor a 3
PBC=0;%condiciones periodicas
dim_site=(N_max+1);%dimension de un sitio


dim=dim_site.^sites;% dimension de todo el espacio
A_id=cell(1,sites);

for i=(0:sites-1)
    A_id{i+1}=eye(dim_site.^i);
end


A_dag=diag(sqrt(1:N_max),+1);
A=A_dag';
N=diag(sqrt(0:N_max));

disp(A_dag)
disp(A)
disp(N)

Ham=BOSE_HUBBARD(A_dag,A,N,A_id,J,U,mu,dim,sites,PBC);
disp(Ham)