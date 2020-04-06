APP=MultMatriz

.PHONY: build clear init run-naming_service run-server run-client

build: ${APP}.idl ${APP}.mpc
	tao_idl -Sa -St ${APP}.idl
	mwc.pl -type gnuace
	make
	@echo "¡Proyecto generado correctamente!"

init:
	tao_idl -Sa -St -GI ${APP}.idl
	@echo "¡Proyecto inicializado correctamente!"

clear:
	rm -f -R *C.* *S.* *Server *Client GNUmakefile* .depend.* .obj/

run-naming_service:
	$TAO_ROOT/orbsvcs/Naming_Service/tao_cosnaming -m 1

run-server: build
	./${APP}Server

run-client: build
	./${APP}Client
