CFLAGS          =
FFLAGS          =
CPPFLAGS        =
FPPFLAGS        =
LOCDIR          = 
EXAMPLESF       = ex1.c ex1f.F90 ex2a.c ex2b.c ex2b_MMS.c ex2b_validation.c ex2e.c
MANSEC          = TS
DIRS            = toy-problem
CLEANFILES      = SA-data/*
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

all: ex1 ex1f ex2a ex2b ex2b_MMS ex2b_validation ex2e
.PHONY: all
