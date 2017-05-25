module typedef
    type surb!surb=surface+subdivision
        real,allocatable::point(:,:)  
        integer,allocatable::line(:,:)
        integer,allocatable::area(:,:)
    end type
    type(surb)::s  !读取初始信息
    real::Vf(20000,3),Ve(20000,3),Vv(20000,3)  
    !vf(新面点)，ve(新边点)，vv(新顶点)
    integer::Lef(30000,2),Lev(30000,2),Na(30000,4) 
    !lef(新面点和新边点连接)，lev(新变点和新顶点连接)，na(新面)
    integer::points,lines,areas,nls
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
    call get_edgepoint()
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
    allocate(s%area(areas,4))
    do i=1,points
        read(10,*)nump,s%point(i,1),s%point(i,2),s%point(i,3) 
    end do
    do j=1,lines
        read(10,"(I5,9X,2(3X,I5))")numl,s%line(j,1),s%line(j,2)
    end do
    do k=1,areas
        read(10,"(I5,9X,4(3X,I5))")numa,s%area(k,1),s%area(k,2),s%area(k,3),s%area(k,4)
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
        do j=1,4
            V=s%point(s%area(i,j),:)
            aveV=aveV+V/4.0           
        end do
        Vf(i,:)=aveV
        V=0.0
        aveV=0.0
    end do
    return
    end subroutine  

!计算新边点
subroutine get_edgepoint()
    use typedef
    implicit none
    real::V(3)
    integer::i,j,k,l,count=0
    do i=1,lines  !线循环
      do j=1,areas
          do k=1,4
              if(s%area(j,k)==s%line(i,1))then
                  do l=1,4
                      if(s%area(j,l)==s%line(i,2).and.k/=l)then
                          v=v+vf(j,:)
                          count=count+1
                      end if
                  end do
              end if
          end do
      end do
      ve(i,:)=(s%point(s%line(i,1),:)+s%point(s%line(i,2),:)+v)/(2+count)
      count=0
      V=0.0
    end do
    return
    end subroutine
    
!计算新顶点
subroutine get_vertexpoint()
    use typedef
    implicit none
    real::Vl(3)=0.0,Va(3)=0.0
    integer::i,j,k,l,m
    real::n=0.0    !n为计数
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
        do l=1,areas
            do m=1,4
                if(s%area(l,m)==i)then
                    Va=Va+Vf(l,:)
                end if
            end do
        end do
        Vv(i,:)=((n-2.0)/n)*s%point(i,:)+Vl/(n**2)+Va/(n**2)
        Vl=0.0
        Va=0.0
        n=0.0
    end do
    return
    end subroutine

!连接
subroutine get_line()
    use typedef
    implicit none
    integer::sx1(4)=(/1,2,3,4/),sx2(4)=(/2,3,4,1/)
    integer::i,j,k,l,m,n,c
    nls=0
    do i=1,areas !新面点循环   Vf——Ve
        do l=1,4
            m=s%area(i,sx1(l))
            n=s%area(i,sx2(l))
          do j=1,lines
            do k=1,2
            if(s%line(j,k)==m)then
                if(k==1)then
                    if(s%line(j,2)==n)then
                       lef(((i-1)*4+l),1)=i
                       lef(((i-1)*4+l),2)=areas+j
                    end if
                else
                    if(s%line(j,1)==n)then
                       lef(((i-1)*4+l),1)=i
                       lef(((i-1)*4+l),2)=areas+j
                    end if
                end if
            end if
            end do
          end do         
       end do
    end do
    
    do k=1,points   !新顶点循环  Vv——Ve
        do l=1,lines
            do m=1,2
                if(s%line(l,m)==k)then
                    nls=nls+1
                    Lev(nls,:)=(/k+areas+lines,l+areas/) !Lev(新顶点，新边点)
                end if
            end do
        end do
    end do
    return
    end subroutine
    
!找出新面的边
subroutine get_newface   
    use typedef
    implicit none
    integer::sx1(5)=(/1,2,3,4,1/),sx2(5)=(/2,3,4,1,2/)
    integer::i,j,k,l,m,n,f
    do i=1,areas !新面点循环   
        do l=1,5
            m=s%area(i,sx1(l))
            n=s%area(i,sx2(l))
          do j=1,lines
            do k=1,2
            if(s%line(j,k)==m)then
                if(k==1)then
                    if(s%line(j,2)==n)then
                        if(l==1)then
                       Na(((i-1)*4+l),1)=i
                       Na(((i-1)*4+l),2)=areas+j
                       Na(((i-1)*4+l),3)=s%area(i,sx2(l))+areas+lines
                        else if(l==5)then
                       Na(((i-1)*4+l-1),4)=areas+j
                        else
                       Na(((i-1)*4+l-1),4)=areas+j 
                       Na(((i-1)*4+l),1)=i
                       Na(((i-1)*4+l),2)=areas+j
                       Na(((i-1)*4+l),3)=s%area(i,sx2(l))+areas+lines
                        end if
                    end if
                else
                    if(s%line(j,1)==n)then
                        if(l==1)then
                       Na(((i-1)*4+l),1)=i
                       Na(((i-1)*4+l),2)=areas+j
                       Na(((i-1)*4+l),3)=s%area(i,sx2(l))+areas+lines
                        else if(l==5)then
                       Na(((i-1)*4+l-1),4)=areas+j
                        else
                       Na(((i-1)*4+l-1),4)=areas+j 
                       Na(((i-1)*4+l),1)=i
                       Na(((i-1)*4+l),2)=areas+j
                       Na(((i-1)*4+l),3)=s%area(i,sx2(l))+areas+lines
                        end if
                    end if
                end if
            end if
            end do
          end do 
       end do
    end do
    return
    end subroutine 
    
!输出可视化文件
subroutine output_data()
    use typedef
    implicit none
    integer::i,j,k,l,m,n,f
    write(11,"(A1,3(3X,I5))")"#",areas+lines+points,4*areas+nls,4*areas
    write(11,"(5(2x,I5))")areas+lines+points,8*areas+nls,1,0,0
    do i=1,areas
        write(11,100) i,Vf(i,1),vf(i,2),vf(i,3) 
    end do 
    do j=1,lines
        write(11,100) areas+j,Ve(j,1),Ve(j,2),Ve(j,3)
    end do
    do k=1,points
        write(11,100) areas+lines+k,Vv(k,1),Vv(k,2),Vv(k,3)
    end do
100 format(I5,2X,3(3X,F9.4))
    
    do l=1,areas
        do m=1,4
            write(11,200)(l-1)*4+m, lef(((l-1)*4+m),1),lef(((l-1)*4+m),2)
        end do
    end do
    do n=1,nls
        write(11,200) 4*areas+n, lev(n,1),lev(n,2)
    end do
200 format(I5,3X,"1",1X,"line",2(3X,I5))
    
    do f=1,4*areas
        write(11,300) f, Na(f,1),Na(f,2),Na(f,3),Na(f,4)
    end do
300 format(I5,3X,"1",1x,"quad",4(3X,I5))
    write(11,"(2(8x,I1))")1,1
    write(11,"(A)")"temputure,"
    deallocate(s%point)
    deallocate(s%line)
    deallocate(s%area)
    return
    end subroutine