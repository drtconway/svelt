FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt upgrade -y
RUN apt install -y \
    python3-full \
    python3-pip
RUN pip install sniffles --break-system-packages
