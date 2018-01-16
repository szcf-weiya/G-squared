PKG_CXXFLAGS=$(shell Rscript -e 'Rcpp:::CxxFlags()')
#PKG_LIBS=$(shell Rscript -e 'Rcpp:::LdFlags()')
PKG_LIBS=-L/home/weiya/R/x86_64-pc-linux-gnu-library/library/Rcpp/libs
GSL_LIBS=$(shell gsl-config --libs)
RCPPGSL_LIBS=-L/home/weiya/R/x86_64-pc-linux-gnu-library/library/RcppGSL/libs
RCPPGSL_FLAGS=-I/home/weiya/R/x86_64-pc-linux-gnu-library/library/RcppGSL/include
omp_CFLAGS=-fopenmp

g2_with_extern.so: g2_with_extern.o
	g++ -shared -o $@ $^ \
			$(PKG_LIBS) -lRcpp \
			$(GSL_LIBS) \
			$(RCPPGSL_LIBS) \
			-L/usr/lib/R/lib -lR

g2_with_extern.o: g2_with_extern.cpp
	g++ -I/usr/share/R/include -DNDEBUG \
	    $(PKG_CXXFLAGS) \
			$(RCPPGSL_FLAGS) \
			$(omp_CFLAGS) \
			-fpic -g -O3 -Wall -c $^ -o $@
