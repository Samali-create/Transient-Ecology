     implicit none

     integer :: mm,i,j,i1,j1,j2,j3,ISEED,IRAND,TT,count,nt,tr
	 
     real*8 :: t0,h,k,k11,k12,k13,R10,C10,P10,KK,xc,yc,xp,yp,R0,C0,y1,y2,y3,RAND, &
     R1,C1,P1,k21,k22,k23,k31,k32,k33,k41,k42,k43,R100,C100,P100, &
     transient,xxx,transient_time,z1,z2,z3
     
     ISEED=time()
     CALL SRAND (ISEED)
     
      do i=1,100
        xxx=rand()
      end do
     
     open(123,file='tt_basin_init.dat')
		
     nt=10000000
!      tr= 990000
     
     xc=0.4d0
     yc=2.009d0

     xp=0.08d0
     yp=2.876d0

     R0=0.16129d0
     C0=0.5d0
	 
     KK=1.0d0 
   	 
	do j1=0,50         !!! loop R
	R100 = 0.1d0+0.8d0*rand()
	write(*,*)j1   
	  do j2=0,50       !!! loop C
	  C100 = 0.1d0+0.6d0*rand()    
			
		do j3=0,50     !!! loop P
		P100 = 0.5d0+0.6d0*rand()     

     R10=R100
     C10=C100   
     P10=P100
       
     h=0.01d0
     t0=0.0d0
     
     do mm=1,nt

      k=R10*(1.0d0-(R10/KK))-xc*yc*(C10*R10/(R10+R0))
      k11=k*h
      k=xc*C10*(((yc*R10)/(R10+R0))-1.0d0)-xp*yp*(C10*P10/(C10+C0))
      k12=k*h
      k=xp*P10*( ((yp*C10)/(C10+C0))-1.0d0)
      k13=k*h

      R1=R10+k11/2.0d0
      C1=C10+k12/2.0d0
      P1=P10+k13/2.0d0
      
      k=R1*(1.0d0-(R1/KK))-xc*yc*(C1*R1/(R1+R0))
      k21=k*h
      k=xc*C1*(((yc*R1)/(R1+R0))-1.0d0)-xp*yp*(C1*P1/(C1+C0))
      k22=k*h
      k=xp*P1*( ((yp*C1)/(C1+C0))-1.0d0)
      k23=k*h

      R1=R10+k21/2.0d0
      C1=C10+k22/2.0d0
      P1=P10+k23/2.0d0
      
      k=R1*(1.0d0-(R1/KK))-xc*yc*(C1*R1/(R1+R0))
      k31=k*h
      k=xc*C1*(((yc*R1)/(R1+R0))-1.0d0)-xp*yp*(C1*P1/(C1+C0))
      k32=k*h
      k=xp*P1*( ((yp*C1)/(C1+C0))-1.0d0)
      k33=k*h

      R1=R10+k31
      C1=C10+k32
      P1=P10+k33
      
      k=R1*(1.0d0-(R1/KK))-xc*yc*(C1*R1/(R1+R0))
      k41=k*h
      k=xc*C1*(((yc*R1)/(R1+R0))-1.0d0)-xp*yp*(C1*P1/(C1+C0))
      k42=k*h
      k=xp*P1*( ((yp*C1)/(C1+C0))-1.0d0)
      k43=k*h

      R10=R10+(k11+2.0*k21+2.0*k31+k41)/6.0d0
      C10=C10+(k12+2.0*k22+2.0*k32+k42)/6.0d0
      P10=P10+(k13+2.0*k23+2.0*k33+k43)/6.0d0

      t0=t0+h
      
      if( (P10.gt.(0.0d0-0.001d0)).and.(P10.lt.(0.0d0+0.001d0)) ) then
        if(t0.gt.1000.0d0)then
            transient_time = t0
            !write(123,*) C100,P100,t0 ! transient_time
        end if  
            exit
      end if    
     
	  
      end do ! Ends of loop time iteration mm
	 
      write(123,*)R100, C100, P100, transient_time   ! transient_time
	  
	  
enddo	      !!! loop R
enddo 		  !!! loop C
enddo         !!! loop P

end program
