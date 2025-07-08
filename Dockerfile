FROM rust:1.88-trixie AS build

ENV SSL_CERT_FILE="/usr/local/share/ca-certificates/ZscalerRootCertificate-2048-SHA256.crt"
RUN curl -o $SSL_CERT_FILE http://repo.mcri.edu.au/ZCC/Certs/ZscalerRootCerts/ZscalerRootCertificate-2048-SHA256.crt
RUN update-ca-certificates

ENV HTSLIB_VER="1.22"
ENV HTSLIB_PKG="htslib-$HTSLIB_VER"
RUN curl -o $HTSLIB_PKG.tar.bz2 https://github.com/samtools/htslib/releases/download/$HTSLIB_VER/$HTSLIB_PKG.tar.bz2
RUN tar -xjf $HTSLIB_PKG.tar.bz2
RUN cd $HTSLIB_PKG && ./configure --prefix=/usr/local && make && make install && rm -rf $HTSLIB_PKG $HTSLIB_PKG.tar.bz2
RUN ldd /usr/local/bin/bgzip | grep '=>' | cut -d ' ' -f 3 | tar -chf /needed-files.tar -T -

ADD Cargo.* /source/
ADD src /source/src

RUN cd /source && cargo build --release

FROM debian:trixie

COPY --from=build /source/target/release/svelt /usr/local/bin/
COPY --from=build /usr/local/bin/bgzip /usr/local/bin/
COPY --from=build /usr/local/bin/htsfile /usr/local/bin/
COPY --from=build /usr/local/bin/tabix /usr/local/bin/
COPY --from=build /usr/local/lib/libhts.so.1.22 /usr/local/lib/
COPY --from=build /needed-files.tar /
RUN tar xf needed-files.tar && rm needed-files.tar
