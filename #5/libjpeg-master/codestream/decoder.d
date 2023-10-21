decoder.o decoder.d : decoder.cpp ../codestream/decoder.hpp ../interface/types.hpp \
  ../config.h ../autoconfig.h ../interface/jpgtypes.hpp \
  ../tools/environment.hpp ../interface/tagitem.hpp \
  ../interface/hooks.hpp ../interface/parameters.hpp ../std/stdlib.hpp \
  ../std/setjmp.hpp ../tools/debug.hpp ../std/assert.hpp \
  ../io/bytestream.hpp ../std/string.hpp ../codestream/tables.hpp \
  ../marker/scantypes.hpp ../boxes/box.hpp ../boxes/databox.hpp \
  ../boxes/namespace.hpp ../boxes/parametrictonemappingbox.hpp \
  ../boxes/tonemapperbox.hpp ../boxes/mergingspecbox.hpp \
  ../boxes/superbox.hpp ../boxes/dctbox.hpp ../boxes/alphabox.hpp \
  ../marker/frame.hpp ../codestream/image.hpp
