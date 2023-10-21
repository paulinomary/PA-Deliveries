checksum.o checksum.d : checksum.cpp ../tools/checksum.hpp ../interface/types.hpp \
  ../config.h ../autoconfig.h ../interface/jpgtypes.hpp \
  ../tools/environment.hpp ../interface/tagitem.hpp \
  ../interface/hooks.hpp ../interface/parameters.hpp ../std/stdlib.hpp \
  ../std/setjmp.hpp debug.hpp ../std/assert.hpp
