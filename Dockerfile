FROM danny-wilson/gcat-project
LABEL app="GCAT with LMM library"
LABEL description="General computational analysis tool with linear mixed model library"
LABEL maintainer="Daniel Wilson"
LABEL build-type="From source"
ENV MKDIR /tmp/libgcat_lmm
RUN mkdir $MKDIR
COPY . $MKDIR
WORKDIR $MKDIR
RUN make
RUN mv lib* /usr/lib/
RUN mv src/* /usr/include/gcat/
RUN rm *.o
WORKDIR /home/ubuntu
ENTRYPOINT ["/usr/bin/gcat"]