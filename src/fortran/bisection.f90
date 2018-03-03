! bisection.f90
!
! Bisection solve for function roots
!
!
!

module random_variate

	abstract interface
		function template(x1,x2,x3,x4) result(y)
			real, intent(in) :: x1,x2,x3,x4
			real :: y
		end function template
		function templatem(x1,x2,x3) result(y)
			real, intent(in) :: x1,x2,x3
			real :: y
		end function templatem
	end interface

	real, parameter :: pi=3.14159265
contains

	function fparallel(p,pp,A,u)
		real, intent (in) :: p,pp,A,u
		real :: fparallel
		real :: gammau,gammap, pu
		pu = u/sqrt(1-u**2)
		gammau = sqrt(1+pu**2)
		gammap = sqrt(1+p**2)

		fparallel = (1+A*gammau*gammap)*exp(-A*(p-pu)**2/(gammap*gammau+p*pu+1))
		return
	end function fparallel
	
	function fparallelm(A,u,pp)
		real, intent (in) :: A,u,pp
		real :: pu
		real :: fparallelm
		pu = u/sqrt(1-u**2)
		fparallelm = pu/A*(1+sqrt(u**2+A**2))
		return
	end function fparallelm

	function fdparallel(p,pp,A,u)
		real, intent (in) :: p,pp,A,u
		real :: fdparallel
		
		fdparallel = fparallel(p,pp,A,u)/fparallel(fparallelm(A,u,pp),pp,A,u)-exp(-2.0)
		return
	end function fdparallel

	function foverfp(p,pp,A,u)
		real, intent(in) :: p,pp,A,u
		real :: gamma_p,p_u,gamma_u
		real :: foverfp,f

		gamma_p = sqrt(1+p**2)
		gamma_u = 1/sqrt(1-u**2)
		p_u = u*gamma_u
		f = p/((1+A*gamma_u*gamma_p)*(2*gamma_p))-(2*(gamma_p*gamma_u+p*p_u+1)*A*(p-p_u) &
		-A*(p-p_u)**2*(gamma_u*p/gamma_p+p_u))/(gamma_p*gamma_u+p*p_u+1)**2
		
		foverfp = 1.0/f
		return
	end function foverfp


	function fps(ps,pp,A,u)
		real, intent (in) :: ps,pp,A,u
		real :: fps
		pu = u/sqrt(1-u**2)
		gammau = sqrt(1+pu**2)
		gammap = sqrt(1+pp**2)
		gammas = sqrt(1+ps**2)

		fps = ps*exp(-A*((pp-pu)**2+gammap**2*gammau**2*ps**2)/(gammap*gammau*gammas+pp*pu+1))
		return
	end function fps

	function fpsm(A,u,pp)
		real, intent (in) :: A,u,pp
		real :: fpsm
		real :: gammau,gammap
		gammau = 1/sqrt(1-u**2)
		gammap = sqrt(1+pp**2)
		fpsm = sqrt(1.0/(2*(A*gammau*gammap)**2)*(1+sqrt(1+(2*A*gammau*gammap)**2)))
		return
	end function fpsm
	
	function fdps(p,pp,A,u)
		real, intent (in) :: p,pp,A,u
		real :: fdps
		
		fdps = fps(p,pp,A,u)/fps(fpsm(A,u,pp),pp,A,u)-exp(-2.0)
		return
	end function fdps
	
	function foverfp2(p_s,p_p,A,u)
		real, intent(in) :: p_s,p_p,A,u
		real :: gamma_p,p_u,gamma_u
		real :: foverfp2,f,N,D

		gamma_p = sqrt(1+p_p**2)
		gamma_s = sqrt(1+p_s**2)
		gamma_u = 1/sqrt(1-u**2)
		p_u = u*gamma_u
		
		N = 2*(gamma_u*gamma_p*gamma_s+p_p*p_u+1)*A*gamma_p**2*gamma_u**2*p_s-A*((p_p-p_u)**2 &
		+gamma_p**2*gamma_u**2*p_s**2)*gamma_u*gamma_p*p_s/gamma_s
    	
		D = gamma_u*gamma_p*gamma_s+p_p*p_u+1
    	
		f = 1./p_s - N/D**2
	
		foverfp2 = 1.0/f
		return
	end function foverfp2

	function findpplus(func,funcm,pp,A,u)
		real, intent(in) :: pp,A,u
		real :: findpplus
		real :: pcand
		procedure(template), pointer :: func
		procedure(templatem), pointer :: funcm

		pcand = funcm(A,u,pp)

		do while (func(pcand,pp,A,u) > 0)
			pcand = pcand*2.0
		end do
		findpplus = pcand
		return
	end function findpplus
	
	function findpminus(func,funcm,pp,A,u)
		real, intent(in) :: pp,A,u
		real :: findpminus
		real :: pcand
		procedure(template), pointer :: func
		procedure(templatem), pointer :: funcm

		pcand = -funcm(A,u,pp)

		do while (func(pcand,pp,A,u) > 0)
			pcand = pcand*2.0
		end do
		findpminus = pcand
		return
	end function findpminus


	function bisect(func,f_a,f_b,A,u,pp)
		real, intent(in) :: f_a, f_b, A, u, pp
		real :: left, right
		integer :: flag
		procedure(template), pointer :: func	
	!	write(*,*) "Entered bisect"

		left = f_a
		right = f_b

		bisect = (left+right)*0.5

		do while(abs(func(bisect,pp,A,u)) > 1e-4)
			bisect = (left+right)*0.5
		!	write(*,*) "Entered Loop: ",left,right
			
			if(func(left,pp,A,u) > 0) then
				flag = 1
			else
				flag = 0
			end if
			if(func(bisect,pp,A,u) > 0 .and. flag == 0) then
				right = bisect
			else if(func(bisect,pp,A,u) < 0 .and. flag == 0) then
				left = bisect
			else if(func(bisect,pp,A,u) > 0.and. flag == 1) then
				left = bisect
			else if(func(bisect,pp,A,u) < 0 .and. flag == 1) then
				right = bisect
			end if
			
		end do
		return
	end function bisect

	function generate_variate(func,pm,pplus,pminus,lambdap,lambdam,qp,qm,qmi,pp)
		real, intent(in) :: pm,pplus,pminus,lambdap,lambdam,qp,qm,qmi,pp
		real :: generate_variate
		procedure(template), pointer :: func
		integer :: flag
		flag = 0
		do while(flag == 0)
			call random_number(U)
			call random_number(V)
			if(U <= qm) then
				Y = U/qm
				X = (1-Y)*(pminus+lambdam)+Y*(pplus-lambdap)
				if(V <= (func(X,pp,A,u)/func(pm,pp,A,u))) then
					flag = 1
				end if
			else if (U <= qm+qp) then
				E = -log((U-qm)/qp)
				X = pplus - lambdap*(1-E)
				if(V <= exp(E)*func(X,pp,A,u)/func(pm,pp,A,u)) then
					flag = 1
				end if
			else
				E = -log((U-(qm+qp))/qmi)
				X = pminus+lambdam*(1-E)
				if(V <= exp(E)*func(X,pp,A,u)/func(pm,pp,A,u)) then
					flag = 1
				end if
			end if
		end do
		generate_variate = X
		return
	end function generate_variate
end module random_variate

program test_bisect
	use random_variate
	
	integer :: Npoints
	real :: fa,fb,A,u
	procedure (template), pointer :: func,funcp, func1,func1p
	procedure (templatem), pointer :: funcm,funcmp
	real :: pplus,pminus,lambdap,lambdam
	real :: psplus,psminus,lambdaps,lambdams,pp
	real, dimension(:), allocatable :: p_parallel, p_s, p_s_sin
	integer :: i


	Npoints = 10000
	allocate(p_parallel(Npoints))
	allocate(p_s(Npoints))
	allocate(p_s_sin(Npoints))
	func => fdparallel
	funcp => fdps
	func1 => fparallel
	func1p => fps
	funcm => fparallelm
	funcmp => fpsm
	
	A = 1
	u = 1e-3
	fa = fparallelm(A,u,1.0)
	fb = findpplus(func,funcm,1.0,A,u)
	fbm = findpminus(func,funcm,1.0,A,u)

	pplus = bisect(func,fa,fb,A,u,1.0)
	pminus = bisect(func,fa,fbm,A,u,1.0)
	pm = fparallelm(A,u,pp)

	lambdap = -foverfp(pplus,pp,A,u)
	lambdam = foverfp(pminus,pp,A,u)

	qmi = lambdam/(pplus-pminus)
	qp = lambdap/(pplus-pminus)
	qm = 1 - (qp+qmi)

	write(*,*) "pm=",pm
	write(*,*) "Value of zeros are ",pplus,pminus
	write(*,*) "Values of lambdas are ",lambdap,lambdam
	write(*,*) "Values of qs are ",qp,qmi,qm
	write(*,*) "About to generate variate"

	i = 0
	do while(i < Npoints)	
		if(modulo(i,10) == 0) then
			write(*,*) "Iteration",i
		end if
	
		p_parallel(i+1) = generate_variate(func1,pm,pplus,pminus,lambdap,lambdam,qp,qm,qmi,0.0)
		pp = p_parallel(i+1)
	
		psm = fpsm(A,u,pp)

		faperp = fpsm(A,u,pp)
		fbperp = findpplus(funcp,funcmp,pp,A,u)
		fbmperp = findpminus(funcp,funcmp,pp,A,u)
		psplus = bisect(funcp,faperp,fbperp,A,u,pp)
		psminus = bisect(funcp,faperp,fbmperp,A,u,pp)
		lambdaps = -foverfp2(psplus,pp,A,u)
		lambdams = foverfp2(psminus,pp,A,u)
		
		qmis = lambdams/(psplus-psminus)
		qps = lambdaps/(psplus-psminus)
		qms = 1 - (qps+qmis)

!		write(*,*) "psm=",psm
!		write(*,*) "Value of zeros are ",psplus,psminus
!		write(*,*) "Values of lambdas are ",lambdaps,lambdams
!		write(*,*) "Values of qs are ",qps,qmis,qms
		p_s(i+1) = generate_variate(func1p,psm,psplus,psminus,lambdaps,lambdams,qps,qms,qmis,pp)
		i = i+1
	end do
	p_s = p_s*sqrt(1+p_parallel**2)
	
	! Distribute p_s uniformly over 2\pi
	i = 0
	do while(i<Npoints)
		call random_number(U)
		p_s_sin(i+1) = p_s(i+1)*sin(U*2*pi)
		p_s(i+1) = p_s(i+1)*cos(U*2*pi)
		i = i+1
	end do


	open(unit=1,file="variates_parallel.txt")
	write(1,*) p_parallel
	close(1)
	open(unit=2,file="variates_perp_1.txt")
	write(2,*) p_s
	close(2)
	open(unit=3,file="variates_perp_2.txt")
	write(3,*) p_s_sin
	close(3)
	
	deallocate(p_parallel)
	deallocate(p_s)
	deallocate(p_s_sin)
end program test_bisect
