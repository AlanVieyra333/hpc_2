FROM yongfu/tao_ace as TAO

WORKDIR /app

ADD ./* ./
RUN make

FROM debian:stable-slim

COPY --from=TAO /opt/ACE_wrappers/bin/* /usr/local/bin/
COPY --from=TAO /opt/ACE_wrappers/lib/* /usr/local/lib/
COPY --from=TAO /opt/ACE_wrappers/TAO/orbsvcs/Naming_Service/tao_cosnaming /usr/local/bin/
COPY --from=TAO /app/MultMatrizClient /app/client
COPY --from=TAO /app/MultMatrizServer /app/server

WORKDIR /app

ENV PATH $PATH:/usr/local/bin
ENV LD_LIBRARY_PATH=/usr/local/lib
