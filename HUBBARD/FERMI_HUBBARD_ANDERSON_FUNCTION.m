function [HAMILTONIAN,NUMERO] = FERMI_HUBBARD_ANDERSON_FUNCTION(P,A_id,P_m,J_v,U,mu,dim,N,Np)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Ham=zeros(dim,dim);
Numero=zeros(dim,dim);

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



N_T=N_do+N_up;%Numero de fermiones por sitio


for i=(1:Np)
    ID1=P_m{i};%matriz que tiene dimensiones dim_site^(i-1)
    ID2=A_id{N-i};%matriz que tiene dimension dim_site^(N-i-1)
    
    Ham=Ham-J_v(i)*(kron(C_d_do*P,kron(ID1,kron(C_do,ID2)))... %dimension= dim_site^{1+i-1+1+N-i-1}=dim_site^{N}
            + kron(C_d_up*P,kron(ID1,kron(C_up,ID2)))...
            + kron(P*C_up,kron(ID1,kron(C_d_up,ID2)))...
            + kron(P*C_do,kron(ID1,kron(C_d_do,ID2))) );
        
end

N_up_do=N_up*N_do;

for i=(1:N)
   ID1=A_id{i};
   ID2=A_id{N-i+1};
   Ham=Ham+U*(kron(ID1,kron(N_up_do,ID2)))+mu*kron(ID1,kron(N_T,ID2)) ;
   Numero=Numero+kron(ID1,kron(N_T,ID2));
end

% P_in=P;
% for i=(1:N-3)
% P_in=kron(P_in,P);
% end

%Ham=Ham+PBC*J*(kron(P*(C_do+C_up),kron(P_in,C_d_do+C_d_up))+kron((C_d_do+C_d_up)*P,kron(P_in,C_do+C_up)));
HAMILTONIAN=Ham;
NUMERO=Numero;



end

