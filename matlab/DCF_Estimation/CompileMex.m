function CompileMex(filename)

eval(sprintf('mex -v CFLAGS="\\$CFLAGS -fopenmp" LDFLAGS="\\$LDFLAGS -fopenmp" %s;',filename)) 