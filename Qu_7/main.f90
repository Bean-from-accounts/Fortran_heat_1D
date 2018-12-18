program chaleur
	use parameters
	implicit none
	
	! Introduction de(s) structure(s)
	type(params)::p
	
	! Introduction des variables
	real,dimension(:),allocatable:: T,T_next,X !vec temp (2), pos
	real:: Dt,Dx !Pas de temps, pas d'espace
		
	!instructions
	call ask_param(p) !demande des parametres
	
	allocate(T(p%N+1),T_next(p%N+1),X(p%N+1)) !alloc dyn
	
	Dt=p%T_fin/p%Nt
	Dx=p%L/p%N
	
	call init(p,T) !initialisation de la temperature
	
	call time_stepping(p,T,T_next,X,Dt,Dx)

	call display(p,T,X)
	
	deallocate(T,T_next,X)

end program chaleur
