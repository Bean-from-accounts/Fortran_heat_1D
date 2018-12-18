program chaleur
	use parameters
	implicit none
	
	! Introduction de(s) structure(s)
	type(params)::p
	
	! Introduction des variables
	real,dimension(:),allocatable:: T,T_next,T_a,X !vec temp (2), pos
	real:: Dt,Dx,time !Pas de temps, pas d'espace
	integer:: j
		
	!instructions
	call ask_param(p) !demande des parametres
	
	allocate(T(p%N+1),T_next(p%N+1),T_a(p%N+1),X(p%N+1)) !alloc dyn
	
	Dt=p%T_fin/p%Nt
	Dx=p%L/p%N
	
	call init(p,T) !initialisation de la temperature
	
	if (p%toggle .eqv. .true.) then
		do j=1,p%Nt
			time=j*Dt
			call time_stepping(p,T,T_next,T_a,X,Dt,Dx,time)
			T_next(1)=p%T_1
			T_next(p%N+1)=p%T_init
			T(:)=T_next(:)
			if (mod(j,10)==0) then
				open(11,file='res.dat')
				call display(p,T,T_a,X)
			end if
		end do
		close(11)
	else
		do j=1,p%Nt
			time=j*Dt
			call time_stepping(p,T,T_next,T_a,X,Dt,Dx,time)
			T_next(1)=p%T_1
			T_next(p%N+1)=p%T_init
			T(:)=T_next(:)
		end do
		call display2(p,T,T_a,X)
	end if
	
	deallocate(T,T_next,T_a,X)

end program chaleur
