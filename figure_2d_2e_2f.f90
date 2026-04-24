!This program is for eps_1 vs eps_2 and result: mean transient time (sw)
     implicit none
	 integer,parameter:: n=20, nt=10**7
     integer :: i,j,i1,j1,j2,ISEED,IRAND,TT,count(100),n,nt  
     
real*8 :: t0,h,k,k1(10,1000),R00(100),C00(100),P00(100),R10(100),C10(100),P10(100),KK,&    	   

xc,yc,xp,yp,R0,C0,y1,y2,y3,xx1,eps1,eps2,A1(100,100),sum1(100),sum2(100),m1,p,deg(100),xx,yy,zz,&     

tot_deg,link_no,R1(100),C1(100),P1(100),transient_time(100),xxx,coup(200),avg_TT,prob

      ISEED=time()
      CALL SRAND (ISEED)
	  xc=0.4d0;yc=2.009d0;xp=0.08d0;yp=2.876d0;R0=0.16129d0;C0=0.5d0;KK=1.0d0

      do i1=1,500000
         xx1=rand()
      end do

open(100,file='ps.txt')

  open(10,file='figure_2a.txt')   ! adjacency matrix
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
              
   open(1000,file='initial_condition.dat.dat')   !initial condition
     do i=1,n
       read(1000,*)xx,yy,zz
       R00(i)=xx
       C00(i)=yy
       P00(i)=zz
     end do     
   CLOSE(1000)

   open(90,file='figure_2d_2e_2f.dat')   !epsilon_1 & epsilon_2
     do i=1,73
       read(90,*)xxx
       coup(i)=xxx
     end do     
   CLOSE(90)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do j2=1,1
      print*,j2
       eps1=coup(j2)

         do j1=1,2
        	!write(*,*)j1
           eps2=coup(j1)

     do j=1,n
       R10(j)=R00(j)
       C10(j)=C00(j)
       P10(j)=P00(j)
     end do
                   
     do j=1,n
      count(j)=0
     end do

     h=0.01d0
     t0=0.0d0
     
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

      k=R10(i)*(1.0d0-(R10(i)/KK))-xc*yc*((C10(i)*R10(i))/(R10(i)+R0))
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

      k=R1(i)*(1.0d0-(R1(i)/KK))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0))
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

      k=R1(i)*(1.0d0-(R1(i)/KK))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0))
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
        
      k=R1(i)*(1.0d0-(R1(i)/KK))-xc*yc*((C1(i)*R1(i))/(R1(i)+R0))
      k1(4,3*i-2)=k*h
      k=xc*C1(i)*(((yc*R1(i))/(R1(i)+R0))-1.0d0)-xp*yp*((C1(i)*P1(i))/(C1(i)+C0))+(eps1/deg(i))*sum1(i)
      k1(4,3*i-1)=k*h
      k=xp*P1(i)*(((yp*C1(i))/(C1(i)+C0))-1.0d0)+(eps2/deg(i))*sum2(i)
      k1(4,3*i)=k*h
     end do 
     
!!!total
     do j=1,n
      R10(j)=R10(j)+(k1(1,3*j-2)+2.0d0*k1(2,3*j-2)+2.0d0*k1(3,3*j-2)+k1(4,3*j-2))/6.0d0
      C10(j)=C10(j)+(k1(1,3*j-1)+2.0d0*k1(2,3*j-1)+2.0d0*k1(3,3*j-1)+k1(4,3*j-1))/6.0d0
      P10(j)=P10(j)+(k1(1,3*j)+2.0d0*k1(2,3*j)+2.0d0*k1(3,3*j)+k1(4,3*j))/6.0d0
     end do
     
     t0=t0+h

  do j=1,n
	if(t0.gt.50)then 
      	if((P10(j).gt.(0.0d0-0.001d0)).and.(P10(j).lt.(0.0d0+0.001d0)))then
        	count(j)=count(j)+1

            if(count(j).eq.1)then
          	  transient_time(j)=t0
            end if      
      	end if
	end if
   !end do

   !do j=1,n
   	if((i1.eq.nt).and.(count(j).eq.0))then
    	transient_time(j)=t0
    end if
  end do

end do ! (i1) time loop
		
!!! average transient time:
	 avg_TT=0.0d0
       do i=1,n
        avg_TT=avg_TT+transient_time(i)
       end do
       write(100,*) eps1,eps2, avg_TT/(1.0d0*n)
                    
     end do  !!! j1
 end do  !!! j2

end program
