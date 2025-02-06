program cont
    implicit none
    integer nstep,natoms,mol,type,ok,ibin,total_bin,id
    integer ntype1,ntype2,i,j
    real x,y,z,vx,vy,vz,xcm,ycm,zcm
    real x1,y1,z1,x2,y2,z2
    real dx,dy,dz,dis,radial_bin,radius,cutoff,p
    real, allocatable:: coord(:,:),m(:,:),radial_dist(:,:),contact_prob(:,:)
    character*30 trjfile
    !!!! Simulation specific parameters
    trjfile="trajectory.lammpstrj"
    print*,trjfile
    cutoff=1.2

    !!!!Initilization 
    nstep=0
    open(unit=11,file=trjfile,status='old')
    call read_lammpstrj_header(11,natoms,ok)
    do i = 1,natoms
        read(11,*)id,mol,type,x,y,z
    enddo
    if(ok<0)print*," Error while reading lammpstrj header"
    !!!!! Matrix Allocations !!!
    allocate(coord(3,natoms))
    allocate(m(natoms,natoms))
   allocate(contact_prob(2,natoms))
    m=0;
    close(11)
    open(unit=11,file=trjfile,status='old')
    do
    call read_lammpstrj_header(11,natoms,ok)
    if(ok<0)exit
    nstep=nstep+1
        do i=1,natoms
        read(11,*)id,mol,type,x,y,z
        coord(1,id)=x
        coord(2,id)=y
        coord(3,id)=z
        enddo
    !    if(nstep>150)then
        do i=1,587
            do j=1,587
            x1 = coord(1,i)
            y1 = coord(2,i)
            z1 = coord(3,i)
            x2 = coord(1,j)
            y2 = coord(2,j)
            z2 = coord(3,j)
            dx = x1-x2
            dy = y1-y2
            dz = z1-z2
            dis= sqrt(dx*dx+dy*dy+dz*dz)
            if(dis<=cutoff)m(i,j)=m(i,j)+1  !contact flag is one if in contact
            enddo
        enddo
 !   endif
    enddo
    open(unit=13,file="contact-probability.dat",status="unknown")
     
    do i=1,587
    do j=1,587
    write(13,"(I8,4X,I8,4X,F12.5)")i,j,m(i,j)/real(nstep)
    enddo
    enddo
    endprogram cont

!!!!!!!! Sub routines !!!!!!!
subroutine read_lammpstrj_header(unt,natoms,ok)
        implicit none
        integer :: unt,ok,natoms
        real :: junk,xlo,xhi,ylo,yhi,zlo,zhi
       character::junk_char*4
        read(unt,*,iostat=ok)junk_char
        if(ok<0)return
        read(unt,*)junk
        read(unt,*)junk_char
        read(unt,*)natoms
        read(unt,*)junk_char
        read(unt,*)xlo,xhi
        read(unt,*)ylo,yhi
        read(unt,*)zlo,zhi
        read(unt,*)junk_char
end subroutine read_lammpstrj_header

