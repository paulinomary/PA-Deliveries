environment.o environment.d : environment.cpp environment.hpp ../config.h \
  ../autoconfig.h ../interface/types.hpp ../interface/jpgtypes.hpp \
  ../interface/tagitem.hpp ../interface/hooks.hpp \
  ../interface/parameters.hpp ../std/stdlib.hpp ../std/setjmp.hpp \
  debug.hpp ../std/assert.hpp ../std/stddef.hpp ../std/string.hpp \
  ../tools/debug.hpp
