predictor.o predictor.d : predictor.cpp ../tools/environment.hpp ../config.h \
  ../autoconfig.h ../interface/types.hpp ../interface/jpgtypes.hpp \
  ../interface/tagitem.hpp ../interface/hooks.hpp \
  ../interface/parameters.hpp ../std/stdlib.hpp ../std/setjmp.hpp \
  ../tools/debug.hpp ../std/assert.hpp ../codestream/predictorbase.hpp \
  ../codestream/predictor.hpp
