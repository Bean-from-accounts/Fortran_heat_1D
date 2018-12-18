subroutine ask_param(p)
	use parameters
	implicit none

	!declarations
	type(params),intent(out):: p

	!instructions
	print*, 'donner valeur de D'
	read(*,*) p%D
	print*, 'donner valeur de L'
	read(*,*) p%L
	print*, 'donner valeur de T_fin'
	read(*,*) p%T_fin
	print*, 'donner valeur de T_init'
	read(*,*) p%T_init
	print*, 'donner valeur de T_1'
	read(*,*) p%T_1
	print*, 'donner valeur de N'
	read(*,*) p%N
	print*, 'donner valeur de Nt'
	read(*,*) p%Nt

end subroutine ask_param


subroutine init(p,T)
	use parameters
	implicit none

	!declarations
	type(params),intent(in):: p
	real,dimension(p%N+1),intent(out):: T
	integer:: i
	
	!instructions
	do i=1,p%N+1
		T(i)=p%T_init
	end do

end subroutine init


subroutine time_stepping(p,T,T_next,X,Dt,Dx)
	use parameters
	implicit none
	
	!declarations
	type(params),intent(in):: p
	real,dimension(p%N+1),intent(out)::T,T_next,X
	real,intent(in)::Dt,Dx
	real::R
	integer::i,j
	
	!instructions
	R=p%D*Dt/(Dx**2.)
	do j=1,p%Nt
		T(p%N+1)=p%T_init
		T(1)=p%T_1
		do i=2,p%N
			T_next(i)=R*T(i-1)+(1.-2.*R)*T(i)+R*T(i+1)
		end do
		T_next(1)=p%T_1
		T_next(p%N+1)=p%T_init
		T(:)=T_next(:)
	end do
	
	do i=1,p%N+1
		X(i)=(i-1)*Dx
	end do

end subroutine time_stepping


subroutine display(p,T,X)
	use parameters
	implicit none

	!declarations
	type(params),intent(in):: p
	real,dimension(p%N+1),intent(in)::T,X
	integer::i
	
	!instructions
	do i=1,p%N+1
		print*, 'position', X(i), 'temperature', T(i)
	end do
	
end subroutine display
