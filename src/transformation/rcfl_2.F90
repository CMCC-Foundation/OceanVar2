     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!Subroutine Gaussian Filter!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine rcfl_2(jpj,G,b1,b2,b3,b)
       integer  c,jpj
       real*8   W(jpj)
       real*8   G(jpj)
       real*8   b1,b2,b3,b



     !!!!!!!!!!!!!!!!!!!!!!!
     !!Forward Differences!!  
     !!!!!!!!!!!!!!!!!!!!!!!
      
         W(1) = b * G(1)
         W(2) = b * G(2) + b1 * W(1)
         W(3) = b * G(3) + b1 * W(2) + b2 * W(1)
          do  c = 4,jpj
            W(c) = b * G(c) +b1 * W(c-1) + b2 * W(c-2) + b3 * W(c-3)
         enddo
     

     !!!!!!!!!!!!!!!!!!!!!!!!!
     !!Backward  Differences!!  
     !!!!!!!!!!!!!!!!!!!!!!!!!

       
        G(jpj)   = b * W(jpj)
        G(jpj-1) = b * W(jpj-1) + b1 * G(jpj)
        G(jpj-2) = b * W(jpj-2) + b1 * G(jpj-1) + b2 * G(jpj)
         do  c = jpj-3,1,-1
              G(c) = b * W(c) +  b1 * G(c+1) + b2 * G(c+2) + b3 * G(c+3)
        enddo

      

      end subroutine rcfl_2




