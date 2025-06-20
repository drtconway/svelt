FROM rust:latest

ENV SSL_CERT_FILE="/usr/local/share/ca-certificates/ZscalerRootCertificate-2048-SHA256.crt"
RUN curl -o $SSL_CERT_FILE http://repo.mcri.edu.au/ZCC/Certs/ZscalerRootCerts/ZscalerRootCertificate-2048-SHA256.crt
RUN update-ca-certificates

ENV HTSLIB_VER="1.22"
ENV HTSLIB_PKG="htslib-$HTSLIB_VER"
RUN wget -O $HTSLIB_PKG.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/$HTSLIB_PKG.tar.bz2
RUN tar -xjf $HTSLIB_PKG.tar.bz2
RUN cd $HTSLIB_PKG && ./configure --prefix=/usr/local && make && make install && rm -rf $HTSLIB_PKG $HTSLIB_PKG.tar.bz2

ADD Cargo.* /source/
ADD src /source/src

RUN cd /source && cargo build --release && cp ./target/release/svelt /usr/local/bin && cd / && rm -rf source

