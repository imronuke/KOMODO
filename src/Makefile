# Makefile

FC = gfortran
FFLAGS = -O4 -cpp -Wall -Wextra -std=f2018

OBJ = \
	  mod_data.o \
	  mod_io.o \
	  mod_xsec.o \
	  mod_nodal.o \
	  mod_gpu.o \
	  mod_cmfd.o \
	  mod_th.o \
	  mod_control.o \
	  mod_trans.o \
	  komodo.o

MOD = $(OBJ:o=mod)

.PHONY : default all clean

EXEC = ../komodo

default : all

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) $(MOD)
	$(RM) $(EXEC)

# dependencies list
komodo.o : mod_data.o mod_io.o mod_control.o mod_trans.o
mod_cmfd.o : mod_data.o mod_gpu.o mod_io.o mod_nodal.o
mod_control.o : mod_data.o mod_io.o mod_xsec.o mod_cmfd.o mod_th.o
mod_data.o : 
mod_gpu.o : mod_data.o
mod_io.o : mod_data.o
mod_nodal.o : mod_data.o mod_io.o
mod_th.o : mod_data.o mod_cmfd.o mod_io.o mod_xsec.o
mod_trans.o : mod_data.o mod_io.o mod_xsec.o mod_cmfd.o mod_th.o mod_control.o mod_nodal.o
mod_xsec.o : mod_data.o mod_io.o

