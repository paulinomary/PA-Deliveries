image.o image.d : image.cpp ../tools/environment.hpp ../config.h ../autoconfig.h \
  ../interface/types.hpp ../interface/jpgtypes.hpp \
  ../interface/tagitem.hpp ../interface/hooks.hpp \
  ../interface/parameters.hpp ../std/stdlib.hpp ../std/setjmp.hpp \
  ../tools/debug.hpp ../std/assert.hpp ../tools/checksum.hpp \
  ../marker/scantypes.hpp ../codestream/image.hpp ../marker/frame.hpp \
  ../boxes/databox.hpp ../boxes/box.hpp ../io/bytestream.hpp \
  ../std/string.hpp ../io/memorystream.hpp ../io/checksumadapter.hpp \
  ../codestream/tables.hpp ../boxes/namespace.hpp \
  ../boxes/parametrictonemappingbox.hpp ../boxes/tonemapperbox.hpp \
  ../boxes/mergingspecbox.hpp ../boxes/superbox.hpp ../boxes/dctbox.hpp \
  ../boxes/alphabox.hpp ../codestream/rectanglerequest.hpp \
  ../tools/rectangle.hpp ../control/bitmapctrl.hpp \
  ../interface/imagebitmap.hpp ../std/stddef.hpp \
  ../control/bufferctrl.hpp ../control/blockctrl.hpp \
  ../control/residualbuffer.hpp ../control/blockbitmaprequester.hpp \
  ../control/blockbuffer.hpp ../control/hierarchicalbitmaprequester.hpp \
  ../boxes/checksumbox.hpp
