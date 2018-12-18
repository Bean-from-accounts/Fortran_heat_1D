subroutine ask_param(p)
	use parameters
	implicit none

	!declarations
	type(params),intent(out):: p

	!instructions
	open(10,file='donnees.dat')
	read(10,*) p%D
	read(10,*) p%L
	read(10,*) p%T_fin
	read(10,*) p%T_init
	read(10,*) p%T_1
	read(10,*) p%N
	read(10,*) p%Nt
	read(10,*) p%toggle
	close(10)

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


subroutine time_stepping(p,T,T_next,T_a,X,Dt,Dx,time)
	use parameters
	implicit none
	
	!declarations
	type(params),intent(in):: p
	real,dimension(p%N+1),intent(out)::T,T_next,T_a,X
	real,intent(in)::Dt,Dx,time
	real::R
	integer::i
	
	!instructions
	R=p%D*Dt/(Dx**2.)
	T(p%N+1)=p%T_init
	T(1)=p%T_1
	
	!calcul position + temp analytique
	do i=1,p%N+1
		X(i)=(i-1)*Dx
		!calcul de la temperature analytique
		T_a(i)=p%T_init+(p%T_1-p%T_init)*(1-erf(X(i)/(2*sqrt(p%D*time))))
	end do
	
	!calcul temp numerique
	do i=2,p%N
		T_next(i)=R*T(i-1)+(1.-2.*R)*T(i)+R*T(i+1)
	end do
	
end subroutine time_stepping


subroutine display(p,T,T_a,X)
	use parameters
	implicit none

	!declarations
	type(params),intent(in):: p
	real,dimension(p%N+1),intent(in)::T,T_a,X
	integer::i
	
	!instructions
	do i=1,p%N+1
		write(11,*) X(i), T(i), T_a(i)
	end do
	
	write(11,*) 
	
end subroutine display


subroutine display2(p,T,T_a,X)
	use parameters
	implicit none

	!declarations
	type(params),intent(in):: p
	real,dimension(p%N+1),intent(in)::T,T_a,X
	integer::i
	
	!instructions
	open(11,file='res.dat')
	do i=1,p%N+1	
		write(11,*) X(i), T(i), T_a(i)
	end do
	close(11)
	
end subroutine display2
