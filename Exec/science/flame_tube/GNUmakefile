PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 3

COMP	   = gnu

USE_MPI    = TRUE

USE_GRAV   = TRUE
USE_REACT  = TRUE

USE_ROTATION  = FALSE
USE_DIFFUSION = TRUE

# define the location of the CASTRO top directory
CASTRO_HOME  ?= ../../..

USE_JACOBIAN_CACHING = TRUE
USE_MODEL_PARSER = TRUE
NUM_MODELS := 2

# This sets the EOS directory in $(MICROPHYSICS_HOME)/EOS
EOS_DIR     := helmholtz

# This sets the network directory in $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := aprox13

INTEGRATOR_DIR := VODE

CONDUCTIVITY_DIR := stellar

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
