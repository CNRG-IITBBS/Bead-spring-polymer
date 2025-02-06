program sox9
implicit none
integer::i,j,k,nstep,nframe,n,ok,atmty,atmid,atmmol,natoms,junk_num
real::cbx,cby,cbz,cbx_avg,cby_avg,cbz_avg,xhi,xlo,yhi,ylo,zhi,zlo,x,y,z
character::trj_filename*111,vel_filename*111,trjtype*6,junk_char*4
real cx(3,587),dist(587)

!this section is a trajectory file reader
nframe=0
trj_filename=trim("poly_wt_5.lammpstrj")
print*,trj_filename
open(unit=12,file=trj_filename)
read(12,*)junk_char
read(12,*)junk_num
read(12,*)junk_char
read(12,*)natoms
read(12,*)junk_char
read(12,*)xlo,xhi
read(12,*)ylo,yhi
read(12,*)zlo,zhi
read(12,*)junk_char
print*,'number of atoms is',natoms
close(12)
dist=0.d0
open(unit=12,file=trj_filename)

do 31
    read(12,'(A10)',iostat=ok)junk_char
    if(ok<0)then
            print*,"frames read=",nframe
            exit
    endif
    nframe=nframe+1
    read(12,'(I10)')junk_num
    read(12,*)junk_char
    read(12,*)natoms
    read(12,*)junk_char
    read(12,*)xlo,xhi
    read(12,*)ylo,yhi
    read(12,*)zlo,zhi
    read(12,*)junk_char
     do 32 i=1,natoms
        read(12,*)atmid,atmmol,atmty,x,y,z
        cx(1,atmid) = x
	cx(2,atmid) = y
	cx(3,atmid) = z
     32 enddo
     do i=1,natoms
	cbx = cx(1,i+1) - cx(1,i)
	cby = cx(2,i+1) - cx(2,i)
	cbz = cx(3,i+1) - cx(3,i)
	dist(i) = dist(i)+sqrt(cbx**2 + cby**2 + cbz**2)
     enddo
31 enddo
do i=1,587
     dist(i)=dist(i)/nframe
enddo
open(unit=13, file='dist.dat')
do i = 1, 587
   write(13, "(I8, 4X, F8.3)") i, dist(i)
end do
close(13)

end program sox9
