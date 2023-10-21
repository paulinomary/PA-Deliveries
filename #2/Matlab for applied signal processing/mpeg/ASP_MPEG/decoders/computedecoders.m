clear all

nBits = 16;

load('../codes/imcodes')

dcddecoder = lookupdecoder(dcdcode, nBits, 'verbose');
save imcodes_dcddecoder dcddecoder;
clear dcddecoder;

nzdecoder = lookupdecoder(nzcode, nBits, 'verbose');
save imcodes_nzdecoder nzdecoder;
clear nzdecoder;

vvdecoder = lookupdecoder(vvcode, nBits, 'verbose');
save imcodes_vvdecoder vvdecoder;
clear vvdecoder;

load('../codes/imdiffcodes')

dcddecoder = lookupdecoder(dcdcode, nBits, 'verbose');
save imdiffcodes_dcddecoder dcddecoder;
clear dcddecoder;

nzdecoder = lookupdecoder(nzcode, nBits, 'verbose');
save imdiffcodes_nzdecoder nzdecoder;
clear nzdecoder;

vvdecoder = lookupdecoder(vvcode, nBits, 'verbose');
save imdiffcodes_vvdecoder vvdecoder;
clear vvdecoder;

load('../codes/immotiondiffcodes')

dcddecoder = lookupdecoder(dcdcode, nBits, 'verbose');
save immotiondiffcodes_dcddecoder dcddecoder;
clear dcddecoder;

nzdecoder = lookupdecoder(nzcode, nBits, 'verbose');
save immotiondiffcodes_nzdecoder nzdecoder;
clear nzdecoder;

vvdecoder = lookupdecoder(vvcode, nBits, 'verbose');
save immotiondiffcodes_vvdecoder vvdecoder;
clear vvdecoder;

load('../codes/mvcodes')
mvvdecoder = lookupdecoder(mvvcode, nBits, 'verbose');
save mvvdecoder;
clear mvvdecoder;

mhvdecoder = lookupdecoder(mhvcode, nBits, 'verbose');
save mhvdecoder;
clear mhvdecoder;

load('../codes/mvdcodes')
mdvvdecoder = lookupdecoder(mdvvcode, nBits, 'verbose');
save mdvvdecoder;
clear mdvvdecoder;

mdhvdecoder = lookupdecoder(mdhvcode, nBits, 'verbose');
save mdhvdecoder;
clear mdhvdecoder;

clear all
load('../codes/pcodes')
pdecoder = lookupdecoder(pcode, nBits, 'verbose');
save pdecoder

clear all