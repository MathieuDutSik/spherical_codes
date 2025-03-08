PROGRAM		= CODE_Analysis CODE_Isomorphism

SRCFILES	= POLY_c_cddlib_mpq.cpp
OBJECTS		= $(SRCFILES:.cpp=.o) $(PROGRAM:%=%.o)
LINKOBJECTS	= $(SRCFILES:.cpp=.o)

CC		= g++ -std=c++17 -Wall -Wextra -g -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../polyhedral_common/src_basic -I../polyhedral_common/src_number -I../polyhedral_common/src_matrix -I../polyhedral_common/src_comb -I../polyhedral_common/src_graph -I../polyhedral_common/src_group -I../polyhedral_common/src_poly -I../polyhedral_common/src_dualdesc -I../polyhedral_common/permutalib/src -I.
#CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE) $(GLPK_INCLUDE) $(CDDLIB_INCLUDE) $(MPREAL_INCLUDE)
CFLAGS 		= -I$(GMP_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NAUTY_INCLUDE)
CHASHMAP        = -I../polyhedral_common/basic_common_cpp/sparse-map/include/tsl/ -I../polyhedral_common/basic_common_cpp/robin-map/include/tsl/ -I../polyhedral_common/basic_common_cpp/hopscotch-map/include/tsl/
#CDEFS           = -DTIMINGS -DDEBUG_SPECIFIC -DDEBUG -DDEBUG_MATRIX_GROUP -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DDEBUG_SPECIFIC -DDEBUG -DDEBUG_MATRIX_GROUP -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CDEFS           = -DTIMINGS
#CDEFS           = -DTIMINGS -DDEBUG
CDEFS           =
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)




#LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(MPFR_LINK) $(BOOST_LINK) $(NAUTY_LINK) $(GLPK_LINK) $(CDDLIB_GMP_LINK)
LIBS		= -lm $(GMP_CXX_LINK) -pthread $(BOOST_LINK) $(NAUTY_LINK)


%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CTOTAL) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(LINKOBJECTS) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
