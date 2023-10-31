!> Program to sum up averaged fields computed for statistics and mean field
!! Martin Karp 27/01-23
!> Add functions to calculate terms needed by TKE budgets
!! Shicheng Dai 25/09-23
program postprocess_fluid_stats
  use neko
  use mean_flow
  implicit none
  
  character(len=NEKO_FNAME_LEN) :: inputchar, mesh_fname, stats_fname, mean_fname
  character(len=NEKO_FNAME_LEN) :: file_to_calculate
  type(file_t) :: mean_file, stats_file, output_file, mesh_file
  real(kind=rp) :: start_time
  type(fld_file_data_t) :: stats_data, mean_data
  type(mean_flow_t) :: avg_flow
  type(fluid_stats_t) :: fld_stats
  type(coef_t) :: coef
  type(dofmap_t) :: dof
  type(space_t) :: Xh
  type(mesh_t) :: msh
  type(gs_t) :: gs_h
  type(field_t), pointer :: u, v, w, p
  type(field_t), target :: pp, uu, vv, ww, uv, uw, vw, tmp1, tmp2
  type(field_list_t) :: reynolds, mean_vel_grad
  type(field_list_t) :: derivative1, derivative2, derivative3
  type(field_list_t) :: derivative4, derivative5, derivative6
  type(field_list_t) :: derivative7, derivative8, derivative9
  integer :: argc, i, n, lx
  
  type(field_t), target :: grad1, grad2, grad3
  type(field_t), target :: grad4, grad5, grad6
  type(field_t), target :: grad7, grad8, grad9
  type(field_t), target :: grad10, grad11, grad12
  type(field_t), target :: grad13, grad14, grad15
  type(field_t), target :: grad16, grad17, grad18
!  
  
  argc = command_argument_count()

  if ((argc .lt. 4) .or. (argc .gt. 4)) then
     if (pe_rank .eq. 0) then
        write(*,*) 'Usage: ./postprocess_fluid_stats mesh.nmsh mean_field.fld stats.fld' 
        write(*,*) 'Example command: ./postprocess_fluid_stats box.nmsh mean_field.fld stats.fld derivative1'   
        write(*,*) 'Computes the statstics from the fld files described in mean_field.fld and stats.fld'
        write(*,*) 'Currently we output several fld files'
        write(*,*) 'In mean_vel_grad:'
        write(*,*) 'x-velocity=dudx'
        write(*,*) 'y-velocity=dudy'
        write(*,*) 'z-velocity=dudz'
        write(*,*) 'pressure=dvdx'
        write(*,*) 'temperature=dvdy'
        write(*,*) 's1=dvdz'
        write(*,*) 's2=dwdx'
        write(*,*) 's2=dwdy'
        write(*,*) 's2=dwdz'
        write(*,*) 'The arrangement of other derivative files can be found in readme'
     end if
     stop
  end if
  
  call neko_init 

  call get_command_argument(1, inputchar) 
  read(inputchar, *) mesh_fname
  mesh_file = file_t(trim(mesh_fname))
  call get_command_argument(2, inputchar) 
  read(inputchar, *) mean_fname
  mean_file = file_t(trim(mean_fname))
  call get_command_argument(3, inputchar) 
  read(inputchar, *) stats_fname
  stats_file = file_t(trim(stats_fname))
  call get_command_argument(4, inputchar) 
  read(inputchar, *) file_to_calculate
  
  call mesh_file%read(msh)
   
  call mean_data%init(msh%nelv,msh%offset_el)
  call stats_data%init(msh%nelv,msh%offset_el)
  call mean_file%read(mean_data)
  call stats_file%read(stats_data)
  
  do i = 1,msh%nelv
     lx = mean_data%lx
     msh%elements(i)%e%pts(1)%p%x(1) = mean_data%x%x(linear_index(1,1,1,i,lx,lx,lx))  
     msh%elements(i)%e%pts(2)%p%x(1) = mean_data%x%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(1) = mean_data%x%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(1) = mean_data%x%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(1) = mean_data%x%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(1) = mean_data%x%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(1) = mean_data%x%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(1) = mean_data%x%x(linear_index(lx,lx,lx,i,lx,lx,lx))

     msh%elements(i)%e%pts(1)%p%x(2) = mean_data%y%x(linear_index(1,1,1,i,lx,lx,lx))  
     msh%elements(i)%e%pts(2)%p%x(2) = mean_data%y%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(2) = mean_data%y%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(2) = mean_data%y%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(2) = mean_data%y%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(2) = mean_data%y%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(2) = mean_data%y%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(2) = mean_data%y%x(linear_index(lx,lx,lx,i,lx,lx,lx))

     msh%elements(i)%e%pts(1)%p%x(3) = mean_data%z%x(linear_index(1,1,1,i,lx,lx,lx))  
     msh%elements(i)%e%pts(2)%p%x(3) = mean_data%z%x(linear_index(lx,1,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(3)%p%x(3) = mean_data%z%x(linear_index(1,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(4)%p%x(3) = mean_data%z%x(linear_index(lx,lx,1,i,lx,lx,lx))
     msh%elements(i)%e%pts(5)%p%x(3) = mean_data%z%x(linear_index(1,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(6)%p%x(3) = mean_data%z%x(linear_index(lx,1,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(7)%p%x(3) = mean_data%z%x(linear_index(1,lx,lx,i,lx,lx,lx))
     msh%elements(i)%e%pts(8)%p%x(3) = mean_data%z%x(linear_index(lx,lx,lx,i,lx,lx,lx))
  end do

  call space_init(Xh, GLL, mean_data%lx, mean_data%ly, mean_data%lz)

  dof = dofmap_t(msh, Xh)
  call gs_init(gs_h, dof)
  call coef_init(coef, gs_h)

  call neko_field_registry%add_field(dof, 'u')
  call neko_field_registry%add_field(dof, 'v')
  call neko_field_registry%add_field(dof, 'w')
  call neko_field_registry%add_field(dof, 'p')

  u => neko_field_registry%get_field('u')
  v => neko_field_registry%get_field('v')
  w => neko_field_registry%get_field('w')
  p => neko_field_registry%get_field('p')

  call avg_flow%init(u, v, w, p)
  call fld_stats%init(coef)
  n = mean_data%u%n
  call copy(avg_flow%u%mf%x,mean_data%u%x,n)
  call copy(avg_flow%v%mf%x,mean_data%v%x,n)
  call copy(avg_flow%w%mf%x,mean_data%w%x,n)
  call copy(avg_flow%p%mf%x,mean_data%p%x,n)
  
  
  call copy(fld_stats%stat_fields%fields(1)%f%x,stats_data%p%x,n)
  call copy(fld_stats%stat_fields%fields(2)%f%x,stats_data%u%x,n)
  call copy(fld_stats%stat_fields%fields(3)%f%x,stats_data%v%x,n)
  call copy(fld_stats%stat_fields%fields(4)%f%x,stats_data%w%x,n)
  call copy(fld_stats%stat_fields%fields(5)%f%x,stats_data%t%x,n)
  do i = 6, size(fld_stats%stat_fields%fields)
     call copy(fld_stats%stat_fields%fields(i)%f%x,stats_data%s(i-5)%x,n)
  end do
  


  if (file_to_calculate .eq. 'mean') then
  
      call field_init(uu,dof)
      call field_init(vv,dof)
      call field_init(ww,dof)
      call field_init(uv,dof)
      call field_init(uw,dof)
      call field_init(vw,dof)
      call field_init(pp,dof)
      call field_init(tmp1,dof)
      call field_init(tmp2,dof)
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      
      
      allocate(mean_vel_grad%fields(12))
      mean_vel_grad%fields(1)%f => pp
      mean_vel_grad%fields(2)%f => uu
      mean_vel_grad%fields(3)%f => vv
      mean_vel_grad%fields(4)%f => ww
      mean_vel_grad%fields(5)%f => uv
      mean_vel_grad%fields(6)%f => uw
      mean_vel_grad%fields(7)%f => vw
      mean_vel_grad%fields(8)%f => tmp1
      mean_vel_grad%fields(9)%f => tmp2
      mean_vel_grad%fields(10)%f => grad1
      mean_vel_grad%fields(11)%f => grad2
      mean_vel_grad%fields(12)%f => grad3
      


      
      call fld_stats%post_process(mean_vel_grad=mean_vel_grad)
      !Fix order of gradients
      mean_vel_grad%fields(2)%f => pp
      mean_vel_grad%fields(3)%f => uu
      mean_vel_grad%fields(4)%f => vv
      mean_vel_grad%fields(1)%f => ww
      mean_vel_grad%fields(5)%f => uv
      mean_vel_grad%fields(6)%f => uw
      mean_vel_grad%fields(7)%f => vw
      mean_vel_grad%fields(8)%f => tmp1
      mean_vel_grad%fields(9)%f => tmp2
      mean_vel_grad%fields(10)%f => grad1
      mean_vel_grad%fields(11)%f => grad2
      mean_vel_grad%fields(12)%f => grad3



      if (pe_rank .eq. 0) write(*,*) 'Writing mean velocity gradient into mean_vel_grad'
      output_file = file_t('mean_vel_grad.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(mean_vel_grad, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if
  
  if (file_to_calculate .eq. 'derivative1') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)
      call field_init(grad10,dof)
      call field_init(grad11,dof)
      call field_init(grad12,dof)
  


      allocate(derivative1%fields(12))
      derivative1%fields(1)%f => grad1
      derivative1%fields(2)%f => grad2
      derivative1%fields(3)%f => grad3
      derivative1%fields(4)%f => grad4
      derivative1%fields(5)%f => grad5
      derivative1%fields(6)%f => grad6
      derivative1%fields(7)%f => grad7
      derivative1%fields(8)%f => grad8
      derivative1%fields(9)%f => grad9
      derivative1%fields(10)%f => grad10
      derivative1%fields(11)%f => grad11
      derivative1%fields(12)%f => grad12


      
      call fld_stats%post_process(derivative1=derivative1)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative1'
      output_file = file_t('derivative1_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative1, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if
  
  if (file_to_calculate .eq. 'derivative2') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)


      allocate(derivative2%fields(9))
      derivative2%fields(1)%f => grad1
      derivative2%fields(2)%f => grad2
      derivative2%fields(3)%f => grad3
      derivative2%fields(4)%f => grad4
      derivative2%fields(5)%f => grad5
      derivative2%fields(6)%f => grad6
      derivative2%fields(7)%f => grad7
      derivative2%fields(8)%f => grad8
      derivative2%fields(9)%f => grad9


      
      call fld_stats%post_process(derivative2=derivative2)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative2'
      output_file = file_t('derivative2_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative2, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if
  
  if (file_to_calculate .eq. 'derivative3') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)
      call field_init(grad10,dof)
      call field_init(grad11,dof)
      call field_init(grad12,dof)


      allocate(derivative3%fields(12))
      derivative3%fields(1)%f => grad1
      derivative3%fields(2)%f => grad2
      derivative3%fields(3)%f => grad3
      derivative3%fields(4)%f => grad4
      derivative3%fields(5)%f => grad5
      derivative3%fields(6)%f => grad6
      derivative3%fields(7)%f => grad7
      derivative3%fields(8)%f => grad8
      derivative3%fields(9)%f => grad9
      derivative3%fields(10)%f => grad10
      derivative3%fields(11)%f => grad11
      derivative3%fields(12)%f => grad12


      
      call fld_stats%post_process(derivative3=derivative3)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative3'
      output_file = file_t('derivative3_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative3, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if
  
  if (file_to_calculate .eq. 'derivative4') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)


      allocate(derivative4%fields(9))
      derivative4%fields(1)%f => grad1
      derivative4%fields(2)%f => grad2
      derivative4%fields(3)%f => grad3
      derivative4%fields(4)%f => grad4
      derivative4%fields(5)%f => grad5
      derivative4%fields(6)%f => grad6
      derivative4%fields(7)%f => grad7
      derivative4%fields(8)%f => grad8
      derivative4%fields(9)%f => grad9



      
      call fld_stats%post_process(derivative4=derivative4)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative4'
      output_file = file_t('derivative4_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative4, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if  

  if (file_to_calculate .eq. 'derivative5') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)


      allocate(derivative5%fields(9))
      derivative5%fields(1)%f => grad1
      derivative5%fields(2)%f => grad2
      derivative5%fields(3)%f => grad3
      derivative5%fields(4)%f => grad4
      derivative5%fields(5)%f => grad5
      derivative5%fields(6)%f => grad6
      derivative5%fields(7)%f => grad7
      derivative5%fields(8)%f => grad8
      derivative5%fields(9)%f => grad9

      
      call fld_stats%post_process(derivative5=derivative5)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative5'
      output_file = file_t('derivative5_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative5, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if
  
  if (file_to_calculate .eq. 'derivative6') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)


      allocate(derivative6%fields(9))
      derivative6%fields(1)%f => grad1
      derivative6%fields(2)%f => grad2
      derivative6%fields(3)%f => grad3
      derivative6%fields(4)%f => grad4
      derivative6%fields(5)%f => grad5
      derivative6%fields(6)%f => grad6
      derivative6%fields(7)%f => grad7
      derivative6%fields(8)%f => grad8
      derivative6%fields(9)%f => grad9

      
      call fld_stats%post_process(derivative6=derivative6)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative6'
      output_file = file_t('derivative6_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative6, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if  

  if (file_to_calculate .eq. 'derivative7') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)


      allocate(derivative7%fields(9))
      derivative7%fields(1)%f => grad1
      derivative7%fields(2)%f => grad2
      derivative7%fields(3)%f => grad3
      derivative7%fields(4)%f => grad4
      derivative7%fields(5)%f => grad5
      derivative7%fields(6)%f => grad6
      derivative7%fields(7)%f => grad7
      derivative7%fields(8)%f => grad8
      derivative7%fields(9)%f => grad9

      
      call fld_stats%post_process(derivative7=derivative7)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative7'
      output_file = file_t('derivative7_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative7, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if  

  if (file_to_calculate .eq. 'derivative8') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)


      allocate(derivative8%fields(9))
      derivative8%fields(1)%f => grad1
      derivative8%fields(2)%f => grad2
      derivative8%fields(3)%f => grad3
      derivative8%fields(4)%f => grad4
      derivative8%fields(5)%f => grad5
      derivative8%fields(6)%f => grad6
      derivative8%fields(7)%f => grad7
      derivative8%fields(8)%f => grad8
      derivative8%fields(9)%f => grad9

      
      call fld_stats%post_process(derivative8=derivative8)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative8'
      output_file = file_t('derivative8_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative8, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if 
  
  if (file_to_calculate .eq. 'derivative9') then
  
      call field_init(grad1,dof)
      call field_init(grad2,dof)
      call field_init(grad3,dof)
      call field_init(grad4,dof)
      call field_init(grad5,dof)
      call field_init(grad6,dof)
      call field_init(grad7,dof)
      call field_init(grad8,dof)
      call field_init(grad9,dof)
      call field_init(grad10,dof)
      call field_init(grad11,dof)
      call field_init(grad12,dof)



      allocate(derivative9%fields(12))
      derivative9%fields(1)%f => grad1
      derivative9%fields(2)%f => grad2
      derivative9%fields(3)%f => grad3
      derivative9%fields(4)%f => grad4
      derivative9%fields(5)%f => grad5
      derivative9%fields(6)%f => grad6
      derivative9%fields(7)%f => grad7
      derivative9%fields(8)%f => grad8
      derivative9%fields(9)%f => grad9
      derivative9%fields(10)%f => grad10
      derivative9%fields(11)%f => grad11
      derivative9%fields(12)%f => grad12



      
      call fld_stats%post_process(derivative9=derivative9)



      if (pe_rank .eq. 0) write(*,*) 'Writing gradient into derivative9'
      output_file = file_t('derivative9_.fld')
      write(*,*) 'time = ', stats_data%time
      call output_file%write(derivative9, stats_data%time)
      if (pe_rank .eq. 0) write(*,*) 'Done'
      call neko_finalize
  	
  end if   
  



end program postprocess_fluid_stats
