
dnl Process this file with autoconf to produce a configure script
AC_INIT
AC_CONFIG_SRCDIR([DESCRIPTION])

dnl Check chmod
chmod +x ./inst/workflow/spmd/2d_I/*.sh
chmod +x ./inst/workflow/spmd/2d_noX/*.sh
chmod +x ./inst/workflow/spmd/2d_V/*.sh
chmod +x ./inst/workflow/spmd/3d_I/*.sh
chmod +x ./inst/workflow/spmd/3d_V/*.sh
chmod +x ./inst/workflow/spmd/3d_noX/*.sh
chmod +x ./inst/workflow/simulation/simu_2d/*.sh
chmod +x ./inst/workflow/simulation/simu_2d/script/*.sh
chmod +x ./inst/workflow/timing/*.sh
chmod +x ./vignettes/*.sh

dnl Start to substitute templates
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
