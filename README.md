# tau-lambda index

To generate the locations of all the minimal factors:
./gen_mf [input text path] [output mf path] [$\tau_\ell$] [$\tau_u$] [$\lambda$] [delimiter/terminal symbols]
+ 'delimiter/terminal symbols': those symbols you use to delimit each sample or the terminal symbols (e.g. "#$")

To gengerate the tau-lambda index of the text (need to gen the minimal factors first):
./gen_index [input text path] [mf path] [output index path] [self-index type] [log](optional)
+ 'self-index type':
    - 0: all kinds of the self-indexes
    - 1: r-index
    - 2: LZ77