#!/bin/bash
R -q -e 'devtools::document()' || exit 1 # Update the doc
R -q -e 'devtools::build_manual(path="manual")' || exit 2 # Create the manual
find manual/ -name "myFun_*" -exec mv {} manual/myFun.pdf \; || exit 4 # Replace the manual
R -q -e 'devtools::check(cran=TRUE, manual=TRUE)' || exit 3 # First check using devtools
R CMD build . || exit 5 # Build the package
R CMD check --as-cran `find . -name "myFun_*.tar.gz"` || exit 6 # Second check using R CMD
find . -name "myFun_*.tar.gz" -delete || exit 7 # Remove the tarball
rm -r myFun.Rcheck || exit 8 # Remove the tmp directory
echo "All done!" && exit 0