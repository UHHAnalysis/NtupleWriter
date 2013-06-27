# Package information
LIBRARY = Ntuple
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

# configure FastJet
INCLUDES += -I$(FASTJETDIR)/include
INCLUDES += -I$(FASTJETDIR)/../include

INCLUDES += -I$(SFRAME_DIR)/core
INCLUDES += -I$(SFRAME_DIR)/core/include
INCLUDES += -I$(SFRAME_DIR)/include

INCLUDES += -I$(SFRAME_DIR)/NtupleWriter
INCLUDES += -I$(SFRAME_DIR)/NtupleWriter/include

# Include the generic compilation rules
include $(SFRAME_DIR)/Makefile.common
