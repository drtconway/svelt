FROM rust:latest

ENV SSL_CERT_FILE="/usr/local/share/ca-certificates/ZscalerRootCertificate-2048-SHA256.crt"
RUN curl -o $SSL_CERT_FILE http://repo.mcri.edu.au/ZCC/Certs/ZscalerRootCerts/ZscalerRootCertificate-2048-SHA256.crt
RUN update-ca-certificates

ADD Cargo.* /source/
ADD src /source/src

RUN cd /source && cargo build --release && cp ./target/release/svelt /usr/local/bin && cd / && rm -rf source