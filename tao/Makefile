APP=MultMatriz

.PHONY: build clear run-naming_service run-server run-client

build: ${APP}.idl IDL.mpc
	tao_idl -Sa -St -GI ${APP}.idl
	mwc.pl -type gnuace
	make
	@echo "¡Proyecto inicializado correctamente!"
	
clear:
	rm -f -R *C.* *S.* *I.* *Server *Client GNUmakefile* .depend.* .obj/

run-naming_service:
	$TAO_ROOT/orbsvcs/Naming_Service/tao_cosnaming -m 1

run-server: build
	./${APP}Server

run-client: build
	./${APP}Client
