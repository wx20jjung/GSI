module write_incr
!$$$ module documentation block
!           .      .    .                                       .
! module:   write_incr 
!   prgmmr: Martin       org:                 date: 2019-09-04
!
! abstract: This module contains routines which write out  
!           the atmospheric increment rather than analysis
!
! program history log:
!   2019-09-04 Martin    Initial version.  Based on ncepnems_io
!
! Subroutines Included:
!   sub write_fv3_increment - writes netCDF increment for FV3 global model
!
!$$$ end documentation block

  implicit none
  private
  public write_fv3_increment

  interface write_fv3_increment
     module procedure write_fv3_inc_
  end interface

contains

  subroutine write_fv3_inc_ (grd,sp_a,filename,mype_out,gfs_bundle,sval,ibin)

!$$$  subprogram documentation block
!                .      .    .
! subprogram:    write_fv3_increment
!
!   prgmmr: Martin            org:                 date: 2019-09-04
!
! abstract: This routine takes GSI analysis increments and writes
!           to a netCDF file on the analysis resolution for use by FV3-GFS
!
! program history log:
!   2019-09-04  martin  Initial version. Based on write_atm_nemsio 
!
!   input argument list:
!     filename  - file to open and write to
!     mype_out  - mpi task to write output file
!    gfs_bundle - bundle containing fields on subdomains
!     ibin      - time bin
!
!   output argument list:
!
! attributes:
!   language: f90
!   machines: ibm RS/6000 SP; SGI Origin 2000; Compaq HP
!
!$$$ end documentation block

! !USES:
    use netcdf
    use kinds, only: r_kind,i_kind

    use mpimod, only: mpi_rtype
    use mpimod, only: mpi_comm_world
    use mpimod, only: ierror
    use mpimod, only: mype

    use gridmod, only: strip, rlats, rlons

    use general_commvars_mod, only: load_grid
    use general_specmod, only: spec_vars
    use general_sub2grid_mod, only: sub2grid_info

    use gsi_bundlemod, only: gsi_bundle, gsi_bundlegetpointer
    use control_vectors, only: lupp

    use constants, only: one, fv, rad2deg

    use gsi_4dvar, only: nobs_bins

    implicit none

! !INPUT PARAMETERS:

    type(sub2grid_info), intent(in) :: grd
    type(spec_vars),     intent(in) :: sp_a
    character(len=24),   intent(in) :: filename  ! file to open and write to
    integer(i_kind),     intent(in) :: mype_out  ! mpi task to write output file
    type(gsi_bundle),    intent(in) :: gfs_bundle
    type(gsi_bundle),    intent(in) :: sval(nobs_bins)
    integer(i_kind),     intent(in) :: ibin      ! time bin

!-------------------------------------------------------------------------

    character(len=120) :: my_name = 'WRITE_FV3INCR'

    !real(r_kind),pointer,dimension(:,:) :: sub_ps
    real(r_kind),pointer,dimension(:,:,:) :: sub_u,sub_v,sub_tsen
    real(r_kind),pointer,dimension(:,:,:) :: sub_q,sub_qana,sub_oz
    real(r_kind),pointer,dimension(:,:,:) :: sub_ql, sub_qi

    real(r_kind),dimension(grd%lat2,grd%lon2,grd%nsig) :: sub_dzb,sub_dza
    real(r_kind),dimension(grd%lat2,grd%lon2,grd%nsig) :: sub_prsl
    real(r_kind),dimension(grd%lat2,grd%lon2,grd%nsig+1) :: sub_prsi
    real(r_kind),dimension(grd%lat2,grd%lon2,grd%nsig+1,ibin) :: ges_geopi

    real(r_kind),dimension(grd%lat1*grd%lon1)     :: psm
    real(r_kind),dimension(grd%lat2,grd%lon2,grd%nsig):: sub_dp
    real(r_kind),dimension(grd%lat1*grd%lon1,grd%nsig):: tsensm,prslm, usm, vsm
    real(r_kind),dimension(grd%lat1*grd%lon1,grd%nsig):: dpsm, qsm, qsmana, ozsm
    real(r_kind),dimension(grd%lat1*grd%lon1,grd%nsig):: qism, qlsm 
    real(r_kind),dimension(grd%lat1*grd%lon1,grd%nsig):: dzsm
    real(r_kind),dimension(max(grd%iglobal,grd%itotsub)) :: work1,work2
    real(r_kind),dimension(grd%nlon,grd%nlat-2):: grid, gridrev
    real(r_kind),dimension(grd%nlon) :: deglons
    real(r_kind),dimension(grd%nlat-2) :: deglats
    real(r_kind),dimension(grd%nsig) :: levsout

    integer(i_kind) :: mm1, i, j, k
    integer(i_kind) :: iret, istatus 
    integer(i_kind) :: ncid_out, lon_dimid, lat_dimid, lev_dimid, ilev_dimid
    integer(i_kind) :: lonvarid, latvarid, levvarid, pfullvarid, ilevvarid, &
                       hyaivarid, hybivarid, uvarid, vvarid, delpvarid, delzvarid, &
                       tvarid, sphumvarid, liqwatvarid, o3varid, icvarid
    integer(i_kind) :: dimids3(3),nccount(3),ncstart(3)

!*************************************************************************
!   Initialize local variables
    mm1=mype+1

    istatus=0
    ! TODO CRM - what is the correct index for sval? always 1? related to nfldsig?
    call gsi_bundlegetpointer(sval(1),'tsen', sub_tsen,  iret); istatus=istatus+iret
    !call gsi_bundlegetpointer(sval(1),'tv', sub_tv,  iret); istatus=istatus+iret
    call gsi_bundlegetpointer(gfs_bundle,'q',  sub_qana,   iret); istatus=istatus+iret ! need this one for Tv to Tsen
    call gsi_bundlegetpointer(sval(1),'q',  sub_q,   iret); istatus=istatus+iret
    call gsi_bundlegetpointer(sval(1),'ql',  sub_ql,   iret); istatus=istatus+iret
    call gsi_bundlegetpointer(sval(1),'qi',  sub_qi,   iret); istatus=istatus+iret
    call gsi_bundlegetpointer(sval(1),'oz', sub_oz,  iret); istatus=istatus+iret
    call gsi_bundlegetpointer(sval(1),'u', sub_u, iret); istatus=istatus+iret
    call gsi_bundlegetpointer(sval(1),'v', sub_v, iret); istatus=istatus+iret
    if ( istatus /= 0 ) then
       if ( mype == 0 ) then
         write(6,*) 'write_fv3_incr_: ERROR'
         write(6,*) 'Missing some of the required fields'
         write(6,*) 'Aborting ... '
      endif
      call stop2(999)
    end if
    
    ! Single task writes increment to file
    if ( mype == mype_out ) then
      ! create the output netCDF file
      call nccheck_incr(nf90_create(path=trim(filename)//".nc", cmode=ior(nf90_clobber,nf90_64bit_offset), ncid=ncid_out))
      ! create dimensions based on analysis resolution, not guess
      call nccheck_incr(nf90_def_dim(ncid_out, "lon", grd%nlon, lon_dimid))
      call nccheck_incr(nf90_def_dim(ncid_out, "lat", grd%nlat-2, lat_dimid))
      call nccheck_incr(nf90_def_dim(ncid_out, "lev", grd%nsig, lev_dimid))
      call nccheck_incr(nf90_def_dim(ncid_out, "ilev", grd%nsig+1, ilev_dimid))
      dimids3 = (/ lon_dimid, lat_dimid, lev_dimid /)
      ! create variables
      call nccheck_incr(nf90_def_var(ncid_out, "lon", nf90_real, (/lon_dimid/), lonvarid))
      call nccheck_incr(nf90_def_var(ncid_out, "lat", nf90_real, (/lat_dimid/), latvarid))
      call nccheck_incr(nf90_def_var(ncid_out, "lev", nf90_real, (/lev_dimid/), levvarid))
      call nccheck_incr(nf90_def_var(ncid_out, "pfull", nf90_real, (/lev_dimid/), pfullvarid))
      call nccheck_incr(nf90_def_var(ncid_out, "ilev", nf90_real, (/ilev_dimid/), ilevvarid))
      call nccheck_incr(nf90_def_var(ncid_out, "hyai", nf90_real, (/ilev_dimid/), hyaivarid))
      call nccheck_incr(nf90_def_var(ncid_out, "hybi", nf90_real, (/ilev_dimid/), hybivarid))
      call nccheck_incr(nf90_def_var(ncid_out, "u_inc", nf90_real, dimids3, uvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "v_inc", nf90_real, dimids3, vvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "delp_inc", nf90_real, dimids3, delpvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "delz_inc", nf90_real, dimids3, delzvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "T_inc", nf90_real, dimids3, tvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "sphum_inc", nf90_real, dimids3, sphumvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "liq_wat_inc", nf90_real, dimids3, liqwatvarid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "o3mr_inc", nf90_real, dimids3, o3varid)) 
      call nccheck_incr(nf90_def_var(ncid_out, "icmr_inc", nf90_real, dimids3, icvarid)) 
      ! place global attributes to parallel calc_increment output
      call nccheck_incr(nf90_put_att(ncid_out, nf90_global, "source", "GSI"))
      call nccheck_incr(nf90_put_att(ncid_out, nf90_global, "comment", &
                                      "global analysis increment from write_fv3_increment"))
      ! add units to lat/lon because that's what the calc_increment utility has
      call nccheck_incr(nf90_put_att(ncid_out, lonvarid, "units", "degrees_east"))
      call nccheck_incr(nf90_put_att(ncid_out, latvarid, "units", "degrees_north"))
      ! end the netCDF file definition
      call nccheck_incr(nf90_enddef(ncid_out))
    end if

    ! Strip off boundary points from subdomains
    call strip(sub_tsen  ,tsensm  ,grd%nsig)
    call strip(sub_q   ,qsm   ,grd%nsig)
    call strip(sub_ql  ,qlsm  ,grd%nsig)
    call strip(sub_qi  ,qism  ,grd%nsig)
    call strip(sub_qana   ,qsmana  ,grd%nsig)
    call strip(sub_oz  ,ozsm  ,grd%nsig)
    !call strip(sub_dp  ,dpsm  ,grd%nsig)
    !call strip(sub_prsl,prslm ,grd%nsig)
    call strip(sub_u   ,usm   ,grd%nsig)
    call strip(sub_v   ,vsm   ,grd%nsig)
    !if (lupp) call strip(sub_dza ,dzsm  ,grd%nsig)

    nccount = (/ grd%nlon, grd%nlat-2, 1 /)
    if (mype == mype_out) then
       ! latitudes
       do j=1,grd%nlat-2
          deglats(j) = rlats(j)*rad2deg
       end do
       ! write to file
       call nccheck_incr(nf90_put_var(ncid_out, latvarid, sngl(deglats), &
                         start = (/1/), count = (/grd%nlat-2/)))
       ! longitudes
       do i=1,grd%nlon
          deglons(i) = rlons(i)*rad2deg
       end do
       ! write to file
       call nccheck_incr(nf90_put_var(ncid_out, lonvarid, sngl(deglons), &
                         start = (/1/), count = (/grd%nlon/)))
       ! levels
       do k=1,grd%nsig
         levsout(k) = float(k)
       end do
       ! write to file
       call nccheck_incr(nf90_put_var(ncid_out, levvarid, sngl(levsout), &
                         start = (/1/), count = (/grd%nsig/)))
       ! pfull
       ! ilev
       ! hyai
       ! hybi
    end if
    ! u increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(usm(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, uvarid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
    ! v increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(vsm(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, vvarid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
    ! delp increment
    ! delz increment
    ! Temperature Increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(tsensm(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, tvarid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
    ! specific humidity increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(qsm(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, sphumvarid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
    ! liquid water increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(qlsm(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, liqwatvarid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
    ! ozone increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(ozsm(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, o3varid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
    ! ice mixing ratio increment
    ncstart = (/ 1, 1, grd%nsig /)
    do k=1,grd%nsig
       call mpi_gatherv(qism(1,k),grd%ijn(mm1),mpi_rtype,&
            work1,grd%ijn,grd%displs_g,mpi_rtype,&
            mype_out,mpi_comm_world,ierror)
       if (mype == mype_out) then
          call load_grid(work1,gridrev)
          ! GSI is N->S, we want S->N 
          do j=1,grd%nlat-2
             grid(:,j) = gridrev(:,grd%nlat-1-j)
          end do
          ! write to file
          call nccheck_incr(nf90_put_var(ncid_out, icvarid, sngl(grid), &
                            start = ncstart, count = nccount))
       endif
       ncstart(3) = grd%nsig-k ! GSI is sfc->top, FV3 is top->sfc
    end do
   ! cleanup and exit
   if ( mype == mype_out ) then
      call nccheck_incr(nf90_close(ncid_out))
      write(6,*) "FV3 netCDF increment written, file= "//trim(filename)//".nc"
   end if


  end subroutine write_fv3_inc_

  subroutine nccheck_incr(status)
    use netcdf
    integer, intent (in   ) :: status
    if (status /= nf90_noerr) then
      print *, "fv3_increment netCDF error", trim(nf90_strerror(status))
      call stop2(999)
    end if
  end subroutine nccheck_incr

end module write_incr