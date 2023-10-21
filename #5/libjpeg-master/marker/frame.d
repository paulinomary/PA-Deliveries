frame.o frame.d : frame.cpp ../tools/environment.hpp ../config.h ../autoconfig.h \
  ../interface/types.hpp ../interface/jpgtypes.hpp \
  ../interface/tagitem.hpp ../interface/hooks.hpp \
  ../interface/parameters.hpp ../std/stdlib.hpp ../std/setjmp.hpp \
  ../tools/debug.hpp ../std/assert.hpp ../tools/checksum.hpp \
  ../marker/frame.hpp ../marker/scantypes.hpp ../boxes/databox.hpp \
  ../boxes/box.hpp ../marker/component.hpp ../tools/rectangle.hpp \
  ../io/bytestream.hpp ../std/string.hpp ../io/checksumadapter.hpp \
  ../codestream/tables.hpp ../boxes/namespace.hpp \
  ../boxes/parametrictonemappingbox.hpp ../boxes/tonemapperbox.hpp \
  ../boxes/mergingspecbox.hpp ../boxes/superbox.hpp ../boxes/dctbox.hpp \
  ../boxes/alphabox.hpp ../codestream/image.hpp ../marker/scan.hpp \
  ../interface/imagebitmap.hpp ../std/stddef.hpp \
  ../control/bitmapctrl.hpp ../control/bufferctrl.hpp \
  ../control/lineadapter.hpp ../tools/line.hpp \
  ../control/blockbitmaprequester.hpp ../control/blockbuffer.hpp \
  ../control/blockctrl.hpp ../control/linebitmaprequester.hpp \
  ../control/linebuffer.hpp ../control/hierarchicalbitmaprequester.hpp \
  ../control/blocklineadapter.hpp ../control/linelineadapter.hpp \
  ../control/residualblockhelper.hpp ../boxes/checksumbox.hpp \
  ../dct/dct.hpp
