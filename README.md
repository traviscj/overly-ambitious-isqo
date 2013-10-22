overly-ambitious-isqo
=====================

a c++ implementation of the [isqo algorithm](http://www.optimization-online.org/DB_HTML/2013/05/3855.html).

Limitations
-----------
* Currently only the "exact algorithm" is implemented.
* You need to install AMPL and qpOASES.

Instructions
------------
1. Type `autoreconf -fvi` to generate the configure script
2. Type `./configure --with-qpOASES_aw=$QPOASESDIR --with-amplsolver=$AMPLDIR && make` to run the configure script and make

Have a lot of fun!
