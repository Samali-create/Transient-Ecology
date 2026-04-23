!This program is for msf::
     
     program potential
     implicit none
     
	integer,parameter::nt=10000000,tr=9900000,n=20
    integer :: i,j,i1,j1,j2,j3,l1,k4,TT,count(100),k2,k3,k5,k6,nn
     
    real*8 :: t0,h,k1(10,1000),R10(100),C10(100),P10(100),K,xc,yc,xp,yp,R0,C0,y1,y2,y3, &   
    eps1,eps2,A1(100,100),sum1(100),sum2(100),m1,p,deg(100), &
    R1(100),C1(100),P1(100),coup(100),R00(100),C00(100),P00(100), RR10(100),CC10(100), &
    PP10(100),RR1(100),CC1(100),PP1(100), &
    v1(3000000),znorm1(3000000),gsc1(3000000),cum1(3000000),w0,temp,coup(10000)

    open(100,file='msf_g_new.dat')
     
     xc=0.4d0
     yc=2.009d0

     xp=0.08d0
     yp=2.876d0

     R0=0.16129d0
     C0=0.5d0  

     K=1.0d0;	
	 eps1=0.00001d0 

open(10,file='global_network.dat')
  do i=1,n
    read(10,*)(A1(i,j),j=1,n)
  end do
close(10)

 do i=1,n
  	deg(i)=0.0d0
    do j=1,n
    	deg(i)=deg(i)+A1(i,j)
    end do
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   
open(1000,file='long_initial_condition.dat')
do i1=1,n
    read(1000,*) R00(i1),C00(i1),P00(i1)
end do     
close(1000)

open(200,file='eps_value.dat')
do j1=1,81
	read(200,*) coup(j1)
enddo

do j2=1,81
  write(*,*)j2
  eps2=coup(j2)   !0.00001+0.0001*j2
 
	
     do j=1,n
       R10(j)=R00(j)
       C10(j)=C00(j)
       P10(j)=P00(j)
     end do
          
      nn=(3*n)*(3*n+1)    
            
      do i=3*n+1,nn
         v1(i)=0.0d0
      end do
      do i=1,3*n 
				v1((3*n+1)*i)=1.0d0
				cum1(i)=0.0d0
      end do
            
     h=0.01d0
     t0=0.0d0
     
     do 10 i1=1,nt   !!!time

!      if(mod(i1,100000).eq.0)then
!        write(*,*)i1
!      end if

!!!1st
		sum1=0.0d0;sum2=0.0d0;
    do i=1,n
			do j=1,n
      	sum1(i)=sum1(i)+A1(i,j)*(C10(j)-C10(i))
        sum2(i)=sum2(i)+A1(i,j)*(P10(j)-P10(i))
      end do

      k1(1,3*i-2)=h*(R10(i)*(1.0d0-(R10(i)/K))-xc*yc*((C10(i)*R10(i))/(R10(i)+R0)))
     
      k1(1,3*i-1)=h*(xc*C10(i)*(((yc*R10(i))/(R10(i)+R0))-1.0d0)-&
      xp*yp*((C10(i)*P10(i))/(C10(i)+C0))+(eps1/deg(i))*sum1(i))
      
      k1(1,3*i)=h*(xp*P10(i)*( ((yp*C10(i))/(C10(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i))
    
      R1(i)=R10(i)+k1(1,3*i-2)/2.0d0
      C1(i)=C10(i)+k1(1,3*i-1)/2.0d0
      P1(i)=P10(i)+k1(1,3*i)/2.0d0
     end do

!!!2nd
		 sum1=0.0d0;sum2=0.0d0;
     do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C1(j)-C1(i))
          sum2(i)=sum2(i)+A1(i,j)*(P1(j)-P1(i))
        end do

      k1(2,3*i-2)=h*(R1(i)*(1.0d0-(R1(i)/K))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0)))
      
      k1(2,3*i-1)=h*(xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-&
      xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i))
      
      k1(2,3*i)=h*(xp*P1(i)*( ((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i))
    
      R1(i)=R10(i)+k1(2,3*i-2)/2.0d0
      C1(i)=C10(i)+k1(2,3*i-1)/2.0d0
      P1(i)=P10(i)+k1(2,3*i)/2.0d0
     end do


!!!3rd
		 sum1=0.0d0;sum2=0.0d0;
     do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C1(j)-C1(i))
          sum2(i)=sum2(i)+A1(i,j)*(P1(j)-P1(i))
        end do

      k1(3,3*i-2)=h*(R1(i)*(1.0d0-(R1(i)/K))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0)))
      
      k1(3,3*i-1)=h*(xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-&
      xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i))
      
      k1(3,3*i)=h*(xp*P1(i)*( ((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i))
    
      R1(i)=R10(i)+k1(3,3*i-2)
      C1(i)=C10(i)+k1(3,3*i-1)
      P1(i)=P10(i)+k1(3,3*i)
     end do
     
!!!!4th
		 sum1=0.0d0;sum2=0.0d0;
		 do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C1(j)-C1(i))
          sum2(i)=sum2(i)+A1(i,j)*(P1(j)-P1(i))
        end do

      	k1(4,3*i-2)=h*(R1(i)*(1.0d0-(R1(i)/K))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0)))
     
      	k1(4,3*i-1)=h*(xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-&
      	xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i))
      
      	k1(4,3*i)=h*(xp*P1(i)*( ((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i))
     end do 
     
     
!!!! total
     do j=1,n
      R10(j)=R10(j)+(k1(1,3*j-2)+2.0d0*k1(2,3*j-2)+2.0d0*k1(3,3*j-2)+k1(4,3*j-2))/6.0d0
      C10(j)=C10(j)+(k1(1,3*j-1)+2.0d0*k1(2,3*j-1)+2.0d0*k1(3,3*j-1)+k1(4,3*j-1))/6.0d0
      P10(j)=P10(j)+(k1(1,3*j)+2.0d0*k1(2,3*j)+2.0d0*k1(3,3*j)+k1(4,3*j))/6.0d0
     end do
     
   if(i1.gt.tr)then
!       write(1,*)t0,(R10(j),C10(j),P10(j),j=1,n)

  do 200 j1=0,3*n-1

     do j=1,n
       RR10(j)=v1((3*j-2)*3*n+1+j1)
       CC10(j)=v1((3*j-1)*3*n+1+j1)
       PP10(j)=v1((3*j)*3*n+1+j1)
     end do
!!!1st      
    sum1=0.0d0;sum2=0.0d0;
    do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(CC10(j)-CC10(i))
          sum2(i)=sum2(i)+A1(i,j)*(PP10(j)-PP10(i))
        end do

      k1(1,3*i-2)=h*(( ( 1.0D0-((2.0D0*R10(i))/K) )-((xc*yc*C10(i)*R0)/(R10(i)+R0)**2.0d0 ) )*RR10(i) &
                        -( (xc*yc*R10(i) )/( R10(i)+R0 ) )*CC10(i))
      
      k1(1,3*i-1)=h*(( (xc*yc*C10(i)*R0 )/( R10(i)+R0 )**2.0d0 )*RR10(i) &
       +( xc*(-1.0d0+((yc*R10(i))/( R10(i)+R0 )))-((xp*yp*P10(i)*C0)/(C10(i)+c0)**2.0D0 ))*CC10(i) &          
           -( (xp*yp*C10(i))/( C10(i)+C0 ))*PP10(i)+(eps1/deg(i))*sum1(i))
      
      k1(1,3*i)=h*(( (xp*yp*P10(i)*C0)/(C10(i)+C0)**2.0d0 )*CC10(i) &
            +( xp*(-1.0d0+((yp*C10(i))/( C10(i)+C0 ))) )*PP10(i)+(eps2/deg(i))*sum2(i))
     
     end do   !!! i
       
     do j=1,n
      RR1(j)=RR10(j)+k1(1,3*j-2)/2.0D0
      CC1(j)=CC10(j)+k1(1,3*j-1)/2.0D0
      PP1(j)=PP10(j)+k1(1,3*j)/2.0D0
     end do
     
!!!2nd     
     sum1=0.0d0;sum2=0.0d0;
     do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(CC1(j)-CC1(i))
          sum2(i)=sum2(i)+A1(i,j)*(PP1(j)-PP1(i))
        end do

      k1(2,3*i-2)=h*(( ( 1.0D0-((2.0D0*R10(i))/K) )-((xc*yc*C10(i)*R0)/(R10(i)+R0)**2.0d0 ) )*RR1(i) &
                        -( (xc*yc*R10(i) )/( R10(i)+R0 ) )*CC1(i))
      
      k1(2,3*i-1)=h*(( (xc*yc*C10(i)*R0 )/( R10(i)+R0 )**2.0d0 )*RR1(i) &
       +( xc*(-1.0d0+((yc*R10(i))/( R10(i)+R0 )))-((xp*yp*P10(i)*C0)/(C10(i)+c0)**2.0D0 ))*CC1(i) &          
           -( (xp*yp*C10(i))/( C10(i)+C0 ))*PP1(i)+(eps1/deg(i))*sum1(i))
      
      k1(2,3*i)=h*( ( (xp*yp*P10(i)*C0)/(C10(i)+C0)**2.0d0 )*CC1(i) &
            +( xp*(-1.0d0+((yp*C10(i))/( C10(i)+C0 ))) )*PP1(i)+(eps2/deg(i))*sum2(i))
    
     end do   !!! i
       
     do j=1,n
      RR1(j)=RR10(j)+k1(2,3*j-2)/2.0D0
      CC1(j)=CC10(j)+k1(2,3*j-1)/2.0D0
      PP1(j)=PP10(j)+k1(2,3*j)/2.0D0
     end do

!!!3rd
     sum1=0.0d0;sum2=0.0d0;
      do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(CC1(j)-CC1(i))
          sum2(i)=sum2(i)+A1(i,j)*(PP1(j)-PP1(i))
        end do

      k1(3,3*i-2)=h*(( ( 1.0D0-((2.0D0*R10(i))/K) )-((xc*yc*C10(i)*R0)/(R10(i)+R0)**2.0d0 ) )*RR1(i) &
                        -( (xc*yc*R10(i) )/( R10(i)+R0 ) )*CC1(i))
      
      k1(3,3*i-1)=h*(( (xc*yc*C10(i)*R0 )/( R10(i)+R0 )**2.0d0 )*RR1(i) &
       +( xc*(-1.0d0+((yc*R10(i))/( R10(i)+R0 )))-((xp*yp*P10(i)*C0)/(C10(i)+c0)**2.0D0 ))*CC1(i) &          
           -( (xp*yp*C10(i))/( C10(i)+C0 ))*PP1(i)+(eps1/deg(i))*sum1(i))
      
      k1(3,3*i)=h*(( (xp*yp*P10(i)*C0)/(C10(i)+C0)**2.0d0 )*CC1(i) &
            +( xp*(-1.0d0+((yp*C10(i))/( C10(i)+C0 ))) )*PP1(i)+(eps2/deg(i))*sum2(i))
     
     end do   !!! i
       
     do j=1,n
      RR1(j)=RR10(j)+k1(3,3*j-2)
      CC1(j)=CC10(j)+k1(3,3*j-1)
      PP1(j)=PP10(j)+k1(3,3*j)
     end do
     
 !!!4th    
     sum1=0.0d0;sum2=0.0d0;
    do i=1,n
        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(CC1(j)-CC1(i))
          sum2(i)=sum2(i)+A1(i,j)*(PP1(j)-PP1(i))
        end do

      k1(4,3*i-2)=h*(( ( 1.0D0-((2.0D0*R10(i))/K) )-((xc*yc*C10(i)*R0)/(R10(i)+R0)**2.0d0 ) )*RR1(i) &
                        -( (xc*yc*R10(i) )/( R10(i)+R0 ) )*CC1(i))
      
      k1(4,3*i-1)=h*(( (xc*yc*C10(i)*R0 )/( R10(i)+R0 )**2.0d0 )*RR1(i) &
       +( xc*(-1.0d0+((yc*R10(i))/( R10(i)+R0 )))-((xp*yp*P10(i)*C0)/(C10(i)+c0)**2.0D0 ))*CC1(i) &          
           -( (xp*yp*C10(i))/( C10(i)+C0 ))*PP1(i)+(eps1/deg(i))*sum1(i))

      
     k1(4,3*i)=h*( ( (xp*yp*P10(i)*C0)/(C10(i)+C0)**2.0d0 )*CC1(i) &
            +( xp*(-1.0d0+((yp*C10(i))/( C10(i)+C0 ))) )*PP1(i)+(eps2/deg(i))*sum2(i))
      
     end do   !!! i
       
    do j=1,n
      RR10(j)=RR10(j)+(k1(1,3*j-2)+2.0d0*k1(2,3*j-2)+2.0d0*k1(3,3*j-2)+k1(4,3*j-2))/6.0d0
      CC10(j)=CC10(j)+(k1(1,3*j-1)+2.0d0*k1(2,3*j-1)+2.0d0*k1(3,3*j-1)+k1(4,3*j-1))/6.0d0
      PP10(j)=PP10(j)+(k1(1,3*j)+2.0d0*k1(2,3*j)+2.0d0*k1(3,3*j)+k1(4,3*j))/6.0d0
 
      v1((3*j-2)*3*n+1+j1)=RR10(j)
      v1((3*j-1)*3*n+1+j1)=CC10(j)
      v1((3*j)*3*n+1+j1)=PP10(j)
    enddo

200   continue    !!!! j1

     
     
  znorm1(1)=0.0d0
	do j1=1,3*n
		znorm1(1)=znorm1(1)+v1(3*n*j1+1)**2.0d0
  end do
  znorm1(1)=sqrt(znorm1(1))
      
	do 40 j1=1,3*n
	  v1(3*n*j1+1)=v1(3*n*j1+1)/znorm1(1)
40  continue

  do 80 j1=2,3*n
      
	do 50 k6=1,(j1-1)
	gsc1(k6)=0.0
	   do 55 l1=1,3*n
	     gsc1(k6)=gsc1(k6)+v1(3*n*l1+j1)*v1(3*n*l1+k6)
55        end do
50      end do

      do 60 k2=1,3*n
	do 65 l1=1,(j1-1)
	  v1(3*n*k2+j1)=v1(3*n*k2+j1)-gsc1(l1)*v1(3*n*k2+l1)
65      end do
60    end do

      znorm1(j1)=0.0d0
	do 70 k3=1,3*n
	  znorm1(j1)=znorm1(j1)+v1(3*n*k3+j1)**2.0d0
70      continue
      znorm1(j1)=sqrt(znorm1(j1))
      
	do 85 k4=1,3*n
	   v1(3*n*k4+j1)=v1(3*n*k4+j1)/znorm1(j1)
85      end do      

80      end do

	do 90 k5=1,3*n
	cum1(k5)=cum1(k5)+log(znorm1(k5))/log(2.0d0)
!	write(*,*)cum1(k5)
90    continue
!!!
      endif

       t0=t0+h

 10    continue ! Time Iteration (i1)

      w0=t0-tr*h
      
      do i=1,3*n
         do j=i+1,3*n
            if(cum1(j).GT.cum1(i))then
               temp=cum1(i)
               cum1(i)=cum1(j)
               cum1(j)=temp
            endif
         enddo
      enddo



!      write(*,*)cum1(1)/w0,cum1(2)/w0,cum1(3)/w0,cum1(4)/w0

      write(100,*) eps2,cum1(1)/w0,cum1(2)/w0,cum1(3)/w0
             
 
     end do  !!! j2

 close(90)
 close(100)
 end program 
