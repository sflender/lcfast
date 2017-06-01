CPP := mpicxx.mpich2
CPPFLAGS := -fopenmp -O3 -g

SRCDIR := ./src
OBJDIR := ./obj
INCLUDES := -I $(SRCDIR) \
			-I ../dtk \
			-I ../../trunk/genericio \
			-I../../software/cfitsio \

OBJECTS := $(OBJDIR)/main.o 

#linking to the libraries
xtract : $(OBJECTS) 
	$(CPP) $(CPPFLAGS) $(OBJECTS) \
    -lgsl -lgslcblas \
	-L../../trunk/hacc.build.datastar/libs \
	-L../dtk/lib \
	-L../../software/cfitsio \
	-lcfitsio \
	-lGenericIO \
	-ldtk \
	-o lcfast

#compilation
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp 
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c $(SRCDIR)/main.cpp -o $(OBJDIR)/main.o

clean:
	rm $(OBJDIR)/*.o lcfast