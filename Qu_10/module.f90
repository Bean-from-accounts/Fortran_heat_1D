module parameters
implicit none

type params
	real:: D,L,T_fin !diff coeff, length, final time
	real:: T_init,T_1 !initial temperature, left bound cond
	integer:: N,Nt !nombre subdivisions spatiales/temporelles
	logical:: toggle !sauvegarde frequentielle ou non
end type params

end module parameters
