!This program is for eps_1 vs eps_2 and result: mean transient time (sw)
     
     implicit none

     integer :: i,j,i1,j1,j2,ISEED,IRAND,TT,count(100),n,nt,nk, realization,tr,rp,q
     
     real*8 :: t0,h,k,k1(10,1000),R10(100),C10(100),P10(100),KK(100),xc,yc,xp,yp,R0,C0,y1,y2,y3,xx1,eps1,eps2, &
	 A1(100,100),sum1(100),sum2(100),m1,p,deg(100),xx,yy,zz,error,tot_deg,link_no,R1(100),C1(100),P1(100), &
	 transient_time(100),xxx,coup(200),avg_TT,prob,e1,SE(10000),BS,BS1,R00(132604),C00(132604),P00(132604)
     
     open(18,file='bs_rd_5.dat')

      ISEED=time()
      CALL SRAND (ISEED)

     xc=0.4d0;yc=2.009d0;xp=0.08d0;yp=2.876d0;R0=0.16129d0;C0=0.5d0;
     n=20;nt=10000000;tr=9000000
       
       do i=1,n
         KK(i)=1.0d0
       end do

      do i1=1,500000
         xx1=rand()
      end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(10,file='random_network.dat')   ! adjacency matrix
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

               
   open(1000,file='basin_initial_condition.dat')   !initial condition
     do i=1,132548
       read(1000,*)xx,yy,zz
       R00(i)=xx
       C00(i)=yy
       P00(i)=zz
     end do     
   CLOSE(1000)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

eps1=0.00001d0

do j1=0,10   
	write(*,*)j1
  eps2=0.01d0+0.0009d0*j1  
		   
		   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Basin stability		   
	realization = 100	   
	do nk=1,realization
       write(*,*)nk
 
    do j = 1,20
    rp=irand()
		
	q=mod(rp,132548)
!    write(*,*)rp,q
	
    R10(j) = R00(q)
	C10(j) = C00(q)
	P10(j) = P00(q)
!	write(*,*) nk, R10(j), C10(j), P10(j)
	end do


     h=0.01d0
     t0=0.0d0
     e1=0.0d0
	 
     do i1=1,nt !!!time loop
			!print*,i1
!!!1st
		do i=1,n
      sum1(i)=0.0d0
      sum2(i)=0.0d0

        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C10(j)-C10(i))
        end do

        do j=1,n
          sum2(i)=sum2(i)+A1(i,j)*(P10(j)-P10(i))
        end do

      k=R10(i)*(1.0d0-(R10(i)/KK(i)))-xc*yc*((C10(i)*R10(i))/(R10(i)+R0))
      k1(1,3*i-2)=k*h
      k=xc*C10(i)*(((yc*R10(i))/(R10(i)+R0))-1.0d0)-xp*yp*((C10(i)*P10(i))/(C10(i)+C0))+(eps1/deg(i))*sum1(i)
      k1(1,3*i-1)=k*h
      k=xp*P10(i)*( ((yp*C10(i))/(C10(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i)
      k1(1,3*i)=k*h
     end do  

     do j=1,n
      R1(j)=R10(j)+k1(1,3*j-2)/2.0d0
      C1(j)=C10(j)+k1(1,3*j-1)/2.0d0
      P1(j)=P10(j)+k1(1,3*j)/2.0d0
     end do

!!!2nd
     do i=1,n

      sum1(i)=0.0d0
      sum2(i)=0.0d0

        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C1(j)-C1(i))
        end do

        do j=1,n
          sum2(i)=sum2(i)+A1(i,j)*(P1(j)-P1(i))
        end do

      k=R1(i)*(1.0d0-(R1(i)/KK(i)))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0))
      k1(2,3*i-2)=k*h
      k=xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i)
      k1(2,3*i-1)=k*h
      k=xp*P1(i)*(((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i)
      k1(2,3*i)=k*h
     end do  

     do j=1,n
      R1(j)=R10(j)+k1(2,3*j-2)/2.0d0
      C1(j)=C10(j)+k1(2,3*j-1)/2.0d0
      P1(j)=P10(j)+k1(2,3*j)/2.0d0
     end do

!!!3rd     
     do i=1,n

      sum1(i)=0.0d0
      sum2(i)=0.0d0

        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C1(j)-C1(i))
        end do

        do j=1,n
          sum2(i)=sum2(i)+A1(i,j)*(P1(j)-P1(i))
        end do

      k=R1(i)*(1.0d0-(R1(i)/KK(i)))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0))
      k1(3,3*i-2)=k*h
      k=xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i)
      k1(3,3*i-1)=k*h
      k=xp*P1(i)*( ((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i)
      k1(3,3*i)=k*h
     end do  

     do j=1,n
      R1(j)=R10(j)+k1(3,3*j-2)
      C1(j)=C10(j)+k1(3,3*j-1)
      P1(j)=P10(j)+k1(3,3*j)
     end do

!!!4th     
     do i=1,n
      sum1(i)=0.0d0
      sum2(i)=0.0d0

        do j=1,n
          sum1(i)=sum1(i)+A1(i,j)*(C1(j)-C1(i))
        end do

        do j=1,n
          sum2(i)=sum2(i)+A1(i,j)*(P1(j)-P1(i))
        end do
        
      k=R1(i)*(1.0d0-(R1(i)/KK(i)))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0))
      k1(4,3*i-2)=k*h
      k=xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i)
      k1(4,3*i-1)=k*h
      k=xp*P1(i)*( ((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i)
      k1(4,3*i)=k*h
     end do 
     
!!!total
     do j=1,n
      R10(j)=R10(j)+(k1(1,3*j-2)+2.0d0*k1(2,3*j-2)+2.0d0*k1(3,3*j-2)+k1(4,3*j-2))/6.0d0
      C10(j)=C10(j)+(k1(1,3*j-1)+2.0d0*k1(2,3*j-1)+2.0d0*k1(3,3*j-1)+k1(4,3*j-1))/6.0d0
      P10(j)=P10(j)+(k1(1,3*j)+2.0d0*k1(2,3*j)+2.0d0*k1(3,3*j)+k1(4,3*j))/6.0d0
     end do
     
     t0=t0+h


		if(i1.gt.tr) then			 	
			 	error=0.0
			 	do i=2,n
			 		error=error+sqrt((R10(i)-R10(1))**2.0+(C10(i)-C10(1))**2.0+(P10(i)-P10(1))**2.0)
			 	enddo
			 		error=error/(n-1)
			 		e1=e1+error
		endif
  end do ! ends (i1)
  e1=e1/((nt-tr)*1.0d0)
		SE(nk)=e1

     end do ! nk                    

     BS=0.0d0
     do nk=1,realization
	    if(SE(nk).lt. 0.001d0)then
		   BS=BS+1.0d0
		end if
     end do	
     BS1= BS/(realization*1.0d0)

     write(18,*)eps2,BS1,BS	 
		
     end do  !!! j1

      end program
