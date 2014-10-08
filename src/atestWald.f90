
!==============================================
!============== Test de Wald multivarie
!==============================================
!b: vecteur de taille m paramètre estimé
!L: matrice indicatrice
!hess :var(b) estime

    subroutine waldmultiv(b,m,L,hess,ddl,wald)

    use optim

    implicit none
    
    integer,intent(in)::m,ddl
    double precision,dimension(m),intent(in)::b
    double precision,dimension(m,m),intent(in)::hess
    integer,dimension(ddl,m),intent(in)::L
    double precision,intent(out)::wald
    double precision,dimension(m,1)::bb
    double precision,dimension(ddl,ddl)::res1
    double precision,dimension(ddl*(ddl+3)/2)::supR
    double precision,dimension(1,1)::resultat
    integer::i,j,ier
    double precision::ep=1.d-30

     bb(:,1)=b
    
    res1 = matmul(L,matmul(hess,transpose(L)))

!stockage de la partie superieure pour utiliser dsinv

    supR = 0.d0
    do i=1,ddl
        do j=i,ddl
            supR((j-1)*j/2+i)=res1(i,j)
        end do
    end do

! inversion de sup_res1

    call dsinvj(supR,ddl,ep,ier)
!    write(*,*)'Critere inversion ',ier

!Reconstruction de la matrice entiere a partir du resultat de dsinv    

    do i=1,ddl
        do j=i,ddl
            res1(i,j)=supR((j-1)*j/2+i)
            res1(j,i)=res1(i,j)
        end do
    end do
    
    do i=2,ddl
        do j=1,i-1
            res1(i,j)=res1(j,i)
        end do
    end do    
    
    resultat = matmul(transpose(matmul(L,bb)),matmul(res1,matmul(L,bb)))

    wald = resultat(1,1)

    end subroutine waldmultiv
