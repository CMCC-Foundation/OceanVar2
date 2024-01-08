   subroutine rcfl_2_ad(jpj,G,b1,b2,b3,b)

       integer       c,jpj
       real*8        W(jpj)
       real*8        G(jpj)
       real*8        b1,b2,b3,b
       real*8        A(jpj)

        W(:)=0
        A(:)=0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Forward  Differences Adjoint !!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do  c =1,jpj-3
           G(c+1)=G(c+1)+b1*G(c)
           G(c+2)=G(c+2)+b2*G(c)
           G(c+3)=G(c+3)+b3*G(c)
           W(c  )=W(c  )+b *G(c)
        enddo

      W(jpj-2)=W(jpj-2)+b*G(jpj-2)
      G(jpj-1)= G(jpj-1)+b1*G(jpj-2)
      G(jpj)=G(jpj)+b2*G(jpj-2)

      W(jpj-1)= W(jpj-1)+b *G(jpj-1)
      G(jpj  )= G(jpj  )+b1*G(jpj-1)

      W(jpj)   = W(jpj)+ b * G(jpj)


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!Backward Differences in Adjoint!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      G(:) = 0.0

        do c =jpj,4,-1
           G(c  )=G(c  )+b *W(c)
           W(c-1)=W(c-1)+b1*W(c)
           W(c-2)=W(c-2)+b2*W(c)
           W(c-3)=W(c-3)+b3*W(c)
        end do

        G(3)=G(3)+b *W(3)
        W(2)=W(2)+b1*W(3)
        W(1)=W(1)+b2*W(3)
  
       G(2)=G(2)+ b*W(2)
       W(1)=W(1)+b1*W(2)
 
       G(1) =G(1)+b * W(1)


    end subroutine rcfl_2_ad
