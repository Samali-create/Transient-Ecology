!This program is for calculating error over global network:

implicit none
integer,parameter::n=20,nt=10000000,tr=9900000
real,parameter::xc=0.4d0,yc=2.009d0,R0=0.16129d0,xp=0.08d0,yp=2.876d0,C0=0.5d0,K=1.0d0

integer :: i,j,i1,j1,mm,a(n,n),i2,j2
real :: R00(n),C00(n),P00(n),R(n),C(n),P(n),R10(n),C10(n),P10(n),&
kR1(n),kR2(n),kR3(n),kR4(n),kC1(n),kC2(n),kC3(n),kC4(n),&
kP1(n),kP2(n),kP3(n),kP4(n),t,h,sum1(n),sum2(n),eps,eps1,eps2,deg(n),&
coup(500),xx,s1,s2,s3,w1,w2,w3,v1,v2,v3,u1,u2,u3,sumr,sumc,sump,error,e1


eps1=0.00001d0;  ! threshold
open(44,file='error.dat')

open(2,file='initial_condition.dat')
 do i=1,n
	read(2,*)R10(i),C10(i),P10(i)
 end do
 close(2)


open(5,file='figure_2a.dat')
 do i=1,n 
	read(5,*) a(i,:)
 end do
 close(5)


do i= 1,n
	deg(i)=0.0d0
	do j=1,n
		deg(i)=deg(i)+a(i,j)
	end do
enddo

open(20,file='figure_2d_2e_2f.dat')
do j1=1,73
   read(20,*)xx
   coup(j1)=xx
 end do
 close(20)
 
 
do i1=1,73 
	eps2=coup(i1)
	print*, i1

 do i=1,n
  R00(i)=R10(i)
  C00(i)=C10(i)
  P00(i)=P10(i)
 enddo

		h=0.01d0
		t=0.0d0
		e1=0.0
		
		do mm=1,nt

!!!!1st
		sum1=0.0d0;sum2=0.0d0
			do i=1,n
				do j=1,n
					sum1(i)=sum1(i)+(a(i,j)*eps1*(C00(j)-C00(i)))/deg(i)
					sum2(i)=sum2(i)+(a(i,j)*eps2*(P00(j)-P00(i)))/deg(i)
				enddo
			
				kR1(i)=R00(i)*(1.0d0-(R00(i)/K))-((xc*yc*C00(i)*R00(i))/(R00(i)+R0))
				kC1(i)=xc*C00(i)*((yc*R00(i))/(R00(i)+R0)-1.0d0)-((xp*yp*C00(i)*P00(i))/(C00(i)+C0))+sum1(i)
				kP1(i)=xp*P00(i)*((yp*C00(i))/(C00(i)+C0)-1.0d0)+sum2(i)
		
				R(i)=R00(i)+h*0.5d0*kR1(i)
				C(i)=C00(i)+h*0.5d0*kC1(i)
				P(i)=P00(i)+h*0.5d0*kP1(i)
			enddo
		
!!!!2nd
		sum1=0.0d0;sum2=0.0d0
			do i=1,n
				do j=1,n
					sum1(i)=sum1(i)+(a(i,j)*eps1*(C(j)-C(i)))/deg(i)
					sum2(i)=sum2(i)+(a(i,j)*eps2*(P(j)-P(i)))/deg(i)
				enddo
			
				kR2(i)=R(i)*(1.0d0-(R(i)/K))-((xc*yc*C(i)*R(i))/(R(i)+R0))
				kC2(i)=xc*C(i)*((yc*R(i))/(R(i)+R0)-1.0d0)-((xp*yp*C(i)*P(i))/(C(i)+C0))+sum1(i)
				kP2(i)=xp*P(i)*((yp*C(i))/(C(i)+C0)-1.0d0)+sum2(i)
		
				R(i)=R00(i)+h*0.5d0*kR2(i)
				C(i)=C00(i)+h*0.5d0*kC2(i)
				P(i)=P00(i)+h*0.5d0*kP2(i)
			enddo
			
!!!!3rd
		sum1=0.0d0;sum2=0.0d0
			do i=1,n
				do j=1,n
					sum1(i)=sum1(i)+(a(i,j)*eps1*(C(j)-C(i)))/deg(i)
					sum2(i)=sum2(i)+(a(i,j)*eps2*(P(j)-P(i)))/deg(i)
				enddo
			
				kR3(i)=R(i)*(1.0d0-(R(i)/K))-((xc*yc*C(i)*R(i))/(R(i)+R0))
				kC3(i)=xc*C(i)*((yc*R(i))/(R(i)+R0)-1.0d0)-((xp*yp*C(i)*P(i))/(C(i)+C0))+sum1(i)
				kP3(i)=xp*P(i)*((yp*C(i))/(C(i)+C0)-1.0d0)+sum2(i)
		
				R(i)=R00(i)+h*kR3(i)
				C(i)=C00(i)+h*kC3(i)
				P(i)=P00(i)+h*kP3(i)
			enddo
		
!!!!	4th
		sum1=0.0d0;sum2=0.0d0
			do i=1,n
				do j=1,n
					sum1(i)=sum1(i)+(a(i,j)*eps1*(C(j)-C(i)))/deg(i)
					sum2(i)=sum2(i)+(a(i,j)*eps2*(P(j)-P(i)))/deg(i)
				enddo
			
				kR4(i)=R(i)*(1.0d0-(R(i)/K))-((xc*yc*C(i)*R(i))/(R(i)+R0))
				kC4(i)=xc*C(i)*((yc*R(i))/(R(i)+R0)-1.0d0)-((xp*yp*C(i)*P(i))/(C(i)+C0))+sum1(i)
				kP4(i)=xp*P(i)*((yp*C(i))/(C(i)+C0)-1.0d0)+sum2(i)
		enddo
			
!!!!	total
			do i=1,n
				R00(i)=R00(i)+(kR1(i)+2.0d0*kR2(i)+2.0d0*kR3(i)+kR4(i))*(h/6.0d0)
				C00(i)=C00(i)+(kC1(i)+2.0d0*kC2(i)+2.0d0*kC3(i)+kC4(i))*(h/6.0d0)
				P00(i)=P00(i)+(kP1(i)+2.0d0*kP2(i)+2.0d0*kP3(i)+kP4(i))*(h/6.0d0)
			enddo
			
			t=t+h	
			
		if(mm.gt.tr) then			 	
			 	error=0.0
			 	do i=2,n
			 		error=error+sqrt((R00(i)-R00(1))**2.0+(C00(i)-C00(1))**2.0+(P00(i)-P00(1))**2.0)
			 	enddo
			 		error=error/(n-1)
			 		e1=e1+error
		endif
			end do   !!! mm
			write(44,*) eps2, e1/(nt-tr)
			
enddo !!!!! i1
 close(44) 

end
