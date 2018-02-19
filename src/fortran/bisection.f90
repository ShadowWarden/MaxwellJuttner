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

	function fps(p,pp,A,u)
		real, intent (in) :: p,pp,A,u
		real :: fperp
		pu = u/sqrt(1-u**2)
		gammau = sqrt(1+pu**2)
		gammap = sqrt(1+pp**2)
		gammas = sqrt(1+p**2)

		fps = ps*exp(-A*((pp-pu)**2+gammap**2*gammau**2*ps**2)/(gammap*gammau*gammas+pp*pu+1))
		return
	end function fps

	function fpsm(A,u,pp)
		real, intent (in) :: A,u,pp
		real :: fpsm
		real :: gammau,gammap
		gammau = 1/sqrt(1-u**2)
		gammap = 1/sqrt(1-pp**2)
		
		fpsm = sqrt(1.0/(2*(A*gammau*gammap)**2)*(1+sqrt(1+(2*A*gammau*gammap)**2)))
		return
	end function fpsm
	
	function fdps(p,pp,A,u)
		real, intent (in) :: p,pp,A,u
		real :: fdps
		
		fdps = fps(p,pp,A,u)/fps(fpsm(A,u,pp),pp,A,u)-exp(-2.0)
		return
	end function fdps

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

end module random_variate

program test_bisect
	use random_variate
	real :: fa,fb,A,u
	procedure (template), pointer :: func
	procedure (templatem), pointer :: funcm

	func => fdparallel
	funcm => fparallelm 
	A = 1
	u = 1e-6
	fa = fparallelm(A,u,1.0)
	fb = findpplus(func,funcm,1.0,A,u)
	fbm = findpminus(func,funcm,1.0,A,u)

	c = bisect(func,fa,fb,A,u,1.0)
	cm = bisect(func,fa,fbm,A,u,1.0)

	write(*,*) "Value of zeros are ",c,cm

end program test_bisect
