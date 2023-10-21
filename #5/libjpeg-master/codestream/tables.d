tables.o tables.d : tables.cpp ../codestream/tables.hpp ../interface/types.hpp \
  ../config.h ../autoconfig.h ../interface/jpgtypes.hpp \
  ../tools/environment.hpp ../interface/tagitem.hpp \
  ../interface/hooks.hpp ../interface/parameters.hpp ../std/stdlib.hpp \
  ../std/setjmp.hpp ../tools/debug.hpp ../std/assert.hpp \
  ../marker/scantypes.hpp ../boxes/box.hpp ../boxes/databox.hpp \
  ../boxes/namespace.hpp ../boxes/parametrictonemappingbox.hpp \
  ../boxes/tonemapperbox.hpp ../boxes/mergingspecbox.hpp \
  ../boxes/superbox.hpp ../boxes/dctbox.hpp ../boxes/alphabox.hpp \
  ../std/string.hpp ../marker/quantization.hpp \
  ../marker/quantizationtable.hpp ../marker/huffmantable.hpp \
  ../marker/component.hpp ../tools/rectangle.hpp ../marker/actable.hpp \
  ../marker/adobemarker.hpp ../marker/restartintervalmarker.hpp \
  ../marker/jfifmarker.hpp ../marker/exifmarker.hpp \
  ../marker/thresholds.hpp ../marker/lscolortrafo.hpp \
  ../marker/frame.hpp ../boxes/inversetonemappingbox.hpp \
  ../boxes/floattonemappingbox.hpp ../boxes/lineartransformationbox.hpp \
  ../boxes/matrixbox.hpp ../boxes/floattransformationbox.hpp \
  ../boxes/checksumbox.hpp ../boxes/filetypebox.hpp \
  ../coding/huffmantemplate.hpp ../coding/actemplate.hpp \
  ../io/bytestream.hpp ../io/checksumadapter.hpp \
  ../colortrafo/colortrafo.hpp ../interface/imagebitmap.hpp \
  ../std/stddef.hpp ../colortrafo/colortransformerfactory.hpp \
  ../tools/traits.hpp ../std/math.hpp ../tools/numerics.hpp \
  ../tools/checksum.hpp ../dct/dct.hpp ../dct/idct.hpp \
  ../dct/liftingdct.hpp
