implicit none
integer,parameter:: nt=2000000, tr=1900000
real:: xc,yc,xp,yp,t, &
kR1,kR2,kR3,kR4,kC1,kC2,kC3,kC4,kP1,kP2,kP3,kP4,&
R0,C0,K,z1,z2,z3,h,C00,R00,P00,R,C,P,fR,fC,fP,x1,x2,x3,y1,y2,y3
integer:: i,j

open (30, file= "p_bif.dat")

do j= 0,1000
	K=0.8+0.0002*j
	!print*, j
	
	R00=0.653;C00=0.427;P00=0.891;h=0.01;t=0.0;

	do i= 1,nt
!!!1st
		kR1=h*fR(R00,C00,P00,K)
		kC1=h*fC(R00,C00,P00)
		kP1=h*fP(R00,C00,P00)
		
		R=R00+(1.0/2.0)*kR1
		C=C00+(1.0/2.0)*kC1
	  	P=P00+(1.0/2.0)*kP1
!!!2nd		
		kR2=h*fR(R,C,P,K)
  		kC2=h*fC(R,C,P)
  		kP2=h*fP(R,C,P)
  	
  		R=R00+(1.0/2.0)*kR2
	  	C=C00+(1.0/2.0)*kC2
	  	P=P00+(1.0/2.0)*kP2
!!!3rd
		kR3=h*fR(R,C,P,K)
	    kC3=h*fC(R,C,P)
	  	kP3=h*fP(R,C,P)
	  
		R=R00+kR3
  		C=C00+kC3
		P=P00+kP3
!!!4th
		kR4=h*fR(R,C,P,K)
	  	kC4=h*fC(R,C,P)
	    kP4=h*fP(R,C,P)
	  
!!!total
		R00=R00+(kR1+2.0*kR2+2.0*kR3+kR4)/6.0
	    C00=C00+(kC1+2.0*kC2+2.0*kC3+kC4)/6.0		
	    P00=P00+(kP1+2.0*kP2+2.0*kP3+kP4)/6.0
	    
		t=t+h
		
		z3 = z2
  		z2 = z1
  		z1 = P00
		
  		
  		if(i.gt.tr)then
  		
			if((z2.gt.z3).and.(z2.gt.z1) )then
				write(30,*) K,z2
				!print*, K,z2
			end if
		endif
		
	end do !!!time
	
	
end do    !!!k loop
 close(30)
end program

real function fR(R,C,P,K)
implicit none
real :: R,C,P,R0,xc,yc,K
xc= 0.4;yc= 2.009;R0= 0.16129!;K=1.0d0;
fR= R*(1.0-(R/K))-((xc*yc*C*R)/(R+R0))
end function

real function fC(R,C,P)
implicit none
real :: R,C,P,R0,C0,xc,yc,xp,yp
xc=0.4;yc= 2.009;xp= 0.08;yp= 2.876;R0=0.16129;C0= 0.5;
fC= xc*C*((yc*R)/(R+R0)-1.0)-((xp*yp*C*P)/(C+C0))
end function

real function fP(R,C,P)
implicit none
real :: R,P,C, C0, xp, yp
xp= 0.08;yp= 2.876;C0= 0.5;
fP= xp*P*((yp*C)/(C+C0)-1.0)
end function
