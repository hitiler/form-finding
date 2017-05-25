module typedef
    type surb!surb=surface+subdivision
        real,allocatable::point(:,:)  
        integer,allocatable::line(:,:)
        integer,allocatable::area(:,:)
    end type
    type(surb)::s  !读取初始信息
    real::Vf(10000,3),Vv(10000,3)  
    !vf(新面点)，vv(新顶点)
    integer::Lfv(10000,2),lff(10000,2),Na(10000,3) 
    !lfv(新面点和新顶点连接)，lff(新顶点和新顶点连接)，na(新面)
    integer::points,lines,areas,lfs,nas   !lfs:新面点连线条数; nas:新面的个数
    end module
    
!**********************************************************
!**********************************************************
!主程序
program xf001
    use typedef
    implicit none
    integer::x
    character(12)::infile,outfile
    do x=1,5
        write(infile,'(a5,i03.3,a4)') "data_",x-1,".inp"
        write(outfile,'(a5,i03.3,a4)') "data_",x,".inp"
    open(10,file=infile,status='old')
    open(11,file=outfile,status='new')
    call input_data()
    call get_facepoint()
    call get_vertexpoint()
    call get_line()
    call get_newface()
    call output_data()
    close(10)
    close(11)
    end do
    stop
    end program xf001    
!**********************************************************  
   
!读取并赋值点，线，面 
subroutine input_data()
    use typedef
    integer::i,j,k,nump,numl,numa
    read(10,"(1x,3(3x,I5))")points,lines,areas !读第一行，点线面的个数
    read(10,*)  !跳过avs信息行
    allocate(s%point(points,3))
    allocate(s%line(lines,2))
    allocate(s%area(areas,3))
    do i=1,points
        read(10,*)nump,s%point(i,1),s%point(i,2),s%point(i,3) 
    end do
    do j=1,lines
        read(10,"(I5,9X,2(3X,I5))")numl,s%line(j,1),s%line(j,2)
    end do
    do k=1,areas
        read(10,"(I5,8X,3(3X,I5))")numa,s%area(k,1),s%area(k,2),s%area(k,3)
    end do    
    return
    end subroutine
  
!计算新面点
subroutine get_facepoint()
    use typedef
    implicit none
    real::V(3),aveV(3)
    integer::i,j
    do i=1,areas  !面循环
        do j=1,3
            V=s%point(s%area(i,j),:)
            aveV=aveV+V/3.0           
        end do
        Vf(i,:)=aveV
!        write(*,*)vf(i,:)
        V=0.0
        aveV=0.0
    end do
    return
    end subroutine  

    
!计算新顶点
subroutine get_vertexpoint()
    use typedef
    implicit none
    real::Vl(3)=0.0
    integer::i,j,k,l,m
    real::an,bn,n=0.0    !n为计数
    real,parameter::pi=3.1415926
    bn=(4.0-2.0*(cos(2.0*pi/2)))/9.0
    do i=1,points  !点循环
        do j=1,lines !循环线来找共边的点
            do k=1,2
                if(s%line(j,k)==i) then !找到了共边的线
                    n=n+1.0                    
                    if(k==1)then
                    Vl=Vl+s%point(s%line(j,2),:)
                    else
                    Vl=Vl+s%point(s%line(j,1),:)
                    end if
                end if
            end do
        end do
        an=(4.0-2.0*(cos(2.0*pi/n)))/9.0
        Vv(i,:)=(1-an)*s%point(i,:)+Vl*an/n
!        write(*,*)n 
!        write(*,*)an
!        write(*,*)vv(i,:)
        Vl=0.0
        n=0.0
    end do
    return
    end subroutine

!连接
subroutine get_line()
    use typedef
    implicit none
    integer::sx1(3)=(/1,2,3/),sx2(3)=(/2,3,1/),sx3(3)=(/3,1,2/),lf(2)
    integer::i,j,k,l,m,n,x,count=0
    lfs=0
    do i=1,areas   !新面点循环  Vf-vv
        do j=1,3
           Lfv((i-1)*3+j,:)=(/i,areas+s%area(i,j)/) !Lfv(新面点，新顶点)
        end do
    end do
    do k=1,lines !新面点循环   Vf——Vf
            m=s%line(k,1)
            n=s%line(k,2)
          do l=1,areas
              do x=1,3
                 if(s%area(l,sx1(x))==m)then
                    if((s%area(l,sx2(x))==n).or.(s%area(l,sx3(x))==n))then
                        count=count+1
                        lf(count)=l
                    end if
                 end if
              end do
          end do
          if(count==1)then
              vv(m,:)=s%point(m,:)
              vv(n,:)=s%point(n,:)
              lfs=lfs+1
              lff(lfs,:)=(/areas+m,areas+n/)
          end if
          if(count==2)then
              lfs=lfs+1
              lff(lfs,:)=lf(:)
          end if
          count=0
          lf=0
    end do
!    write(*,*)lfs
    return
    end subroutine
    
!找出新面的边
subroutine get_newface   
    use typedef
    implicit none
    integer::sx1(3)=(/1,2,3/),sx2(3)=(/2,3,1/),sx3(3)=(/3,1,2/),lf(2)
    integer::i,j,k,l,m,n,x,count=0
    nas=0
    do k=1,lines 
            m=s%line(k,1)
            n=s%line(k,2)
          do l=1,areas
              do x=1,3
                 if(s%area(l,sx1(x))==m)then
                    if((s%area(l,sx2(x))==n).or.(s%area(l,sx3(x))==n))then
                        count=count+1
                        lf(count)=l
                    end if
                 end if
              end do
          end do
          if(count==2)then
              nas=nas+2
              na(nas-1,:)=(/lf(1),m+areas,lf(2)/) 
              na(nas,:)=(/lf(2),n+areas,lf(1)/)!没有判断逆时针
          end if
          if(count==1)then
              nas=nas+1
              na(nas,:)=(/lf(1),m+areas,n+areas/)
          end if
          count=0
    end do
!    write(*,*)nas
    return
    end subroutine 
    
!输出可视化文件
subroutine output_data()
    use typedef
    implicit none
    integer::i,j,k,l,m,n,f
    write(11,"(A1,3(3X,I5))")"#",areas+points,3*areas+lfs,3*areas
    write(11,"(5(2x,I5))")areas+points,6*areas+lfs,1,0,0
    do i=1,areas
        write(11,100) i,Vf(i,1),vf(i,2),vf(i,3) 
    end do
    
    do k=1,points
        write(11,100) areas+k,Vv(k,1),Vv(k,2),Vv(k,3)
    end do
100 format(I5,2X,3(3X,F9.4))
    
    do l=1,areas
        do m=1,3
            write(11,200)(l-1)*3+m, lfv(((l-1)*3+m),1),lfv(((l-1)*3+m),2)
        end do
    end do
    
    do n=1,lfs
        write(11,200) 3*areas+n, lff(n,1),lff(n,2)
    end do
200 format(I5,3X,"1",1X,"line",2(3X,I5))
    
    do f=1,nas
        write(11,300) f, Na(f,1),Na(f,2),Na(f,3)
    end do
300 format(I5,3X,"1",1x,"tri",3(3X,I5))
    write(11,"(2(8x,I1))")1,1
    write(11,"(A)")"temputure,"
    deallocate(s%point)
    deallocate(s%line)
    deallocate(s%area)
    return
    end subroutine