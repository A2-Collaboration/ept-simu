c
      real hist(2,75)
      open(1,file="brel.dat",status="old")
      de=0.06294
      sum=0
      sum1=0.
      sum2=0.
      sum3=0.
      sum4=0.
      do i=1,75
        read(1,*) hist(1,i),hist(2,i)
        sum=sum+hist(2,i)
      end do
      do i=1,75
        hist(2,i)=hist(2,i)/sum
      end do
      do i=1,75
        print*,i,(hist(j,i),j=1,2)
        sum1=sum1+hist(2,i)
        sum2=sum2+hist(1,i)*hist(2,i)
        sum3=sum3+hist(1,i)**2*hist(2,i)
      end do
      print*,'hist(1,50),hist(2,50)',hist(1,50),hist(2,50)
      print*,'sum1,sum2,sum3     ',sum1,sum2,sum3
      print*,'sum2/sum1,sum3/sum1',sum2/sum1,sum3/sum1
      sig=sum3/sum1-(sum2/sum1)**2
      print*,'sig,sqrt(sig)',sig,sqrt(sig)
      do i=1,75
        sum4=sum4+hist(2,i)*(hist(1,i)-sum2/sum1)**2
      end do
      print*,'sum4,sqrt(sum4)',sum4,sqrt(sum4)
      stop
      end
