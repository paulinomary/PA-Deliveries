scan.o scan.d : scan.cpp ../marker/scan.hpp ../tools/environment.hpp ../config.h \
  ../autoconfig.h ../interface/types.hpp ../interface/jpgtypes.hpp \
  ../interface/tagitem.hpp ../interface/hooks.hpp \
  ../interface/parameters.hpp ../std/stdlib.hpp ../std/setjmp.hpp \
  ../tools/debug.hpp ../std/assert.hpp ../tools/rectangle.hpp \
  ../interface/imagebitmap.hpp ../std/stddef.hpp ../marker/scantypes.hpp \
  ../io/bytestream.hpp ../std/string.hpp ../marker/frame.hpp \
  ../boxes/databox.hpp ../boxes/box.hpp ../marker/component.hpp \
  ../codestream/tables.hpp ../boxes/namespace.hpp \
  ../boxes/parametrictonemappingbox.hpp ../boxes/tonemapperbox.hpp \
  ../boxes/mergingspecbox.hpp ../boxes/superbox.hpp ../boxes/dctbox.hpp \
  ../boxes/alphabox.hpp ../codestream/entropyparser.hpp \
  ../codestream/sequentialscan.hpp ../io/bitstream.hpp \
  ../tools/checksum.hpp ../coding/quantizedrow.hpp \
  ../coding/blockrow.hpp ../codestream/acsequentialscan.hpp \
  ../coding/qmcoder.hpp ../codestream/losslessscan.hpp \
  ../codestream/predictivescan.hpp ../codestream/predictorbase.hpp \
  ../tools/line.hpp ../codestream/aclosslessscan.hpp \
  ../codestream/refinementscan.hpp ../io/memorystream.hpp \
  ../codestream/acrefinementscan.hpp \
  ../codestream/singlecomponentlsscan.hpp ../codestream/jpeglsscan.hpp \
  ../control/linebuffer.hpp ../codestream/lineinterleavedlsscan.hpp \
  ../codestream/sampleinterleavedlsscan.hpp \
  ../coding/huffmantemplate.hpp ../marker/huffmantable.hpp \
  ../marker/actable.hpp ../marker/thresholds.hpp \
  ../control/bitmapctrl.hpp ../control/bufferctrl.hpp
