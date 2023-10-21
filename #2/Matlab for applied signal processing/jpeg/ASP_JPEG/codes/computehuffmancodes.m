clear all
load statscodecs3_4
dcdcode = huffmancode(hdcd, 'nowarnings', 'zeroprobcodes');
nzcode = huffmancode(hnz, 'nowarnings', 'zeroprobcodes');
hac(129) = [];
vvcode = huffmancode(hac, 'nowarnings', 'zeroprobcodes');
vvcode = {vvcode{1, 1 : 128}, '!', vvcode{1, 129 : end}};
save imcodes dcdcode nzcode vvcode

clear all
load statscodec5
dcdcode = huffmancode(hdcd, 'nowarnings', 'zeroprobcodes');
nzcode = huffmancode(hnz, 'nowarnings', 'zeroprobcodes');
hac(129) = [];
vvcode = huffmancode(hac, 'nowarnings', 'zeroprobcodes');
vvcode = {vvcode{1, 1 : 128}, '!', vvcode{1, 129 : end}};
save imdiffcodes dcdcode nzcode vvcode

clear all
load statscodec6
dcdcode = huffmancode(hdcd, 'nowarnings', 'zeroprobcodes');
nzcode = huffmancode(hnz, 'nowarnings', 'zeroprobcodes');
hac(129) = [];
vvcode = huffmancode(hac, 'nowarnings', 'zeroprobcodes');
vvcode = {vvcode{1, 1 : 128}, '!', vvcode{1, 129 : end}};
save immotiondiffcodes dcdcode nzcode vvcode

mvvcode = huffmancode(hvv, 'nowarnings', 'zeroprobcodes');
mhvcode = huffmancode(hvh, 'nowarnings', 'zeroprobcodes');
save mvcodes mvvcode mhvcode

mdvvcode = huffmancode(hdvv, 'nowarnings', 'zeroprobcodes');
mdhvcode = huffmancode(hdvh, 'nowarnings', 'zeroprobcodes');
save mvdcodes mdvvcode mdhvcode

clear all
load statspredictivemask
pcode = huffmancode(probsp, 'nowarnings');
save pcode

clear all