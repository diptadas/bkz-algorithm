FROM alpine
RUN apk add --update wget make g++ m4 perl

RUN mkdir /bkz-tools
WORKDIR /bkz-tools

# install GMP
RUN wget https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2 && \
tar xjf gmp-6.1.2.tar.bz2 && cd gmp-6.1.2 && \
./configure && make && make install

# install mpfr
RUN wget https://www.mpfr.org/mpfr-current/mpfr-4.0.1.tar.gz && \
tar xzf mpfr-4.0.1.tar.gz && cd mpfr-4.0.1 && \
./configure && make && make install

# install fplll
RUN wget https://github.com/fplll/fplll/releases/download/5.2.1/fplll-5.2.1.tar.gz && \
tar xzf fplll-5.2.1.tar.gz && cd fplll-5.2.1 && \
./configure && make && make install

# install NTL
RUN wget https://www.shoup.net/ntl/ntl-5.5.2.tar.gz && \
tar xzf ntl-5.5.2.tar.gz && cd ntl-5.5.2/src && \
./configure && make && make install

# cleanup
WORKDIR /
RUN rm -r /bkz-tools