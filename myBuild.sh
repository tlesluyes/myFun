#!/bin/bash
R -q -e 'devtools::document()'
R -q -e 'devtools::build_manual(path="manual")'
find manual/ -name "myFun_*" -exec mv {} manual/myFun.pdf \;
R CMD build .
R CMD check `find . -name "myFun_*.tar.gz"`
find . -name "myFun_*.tar.gz" -delete
rm -r myFun.Rcheck