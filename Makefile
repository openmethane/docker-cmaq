.phony: build
build: Dockerfile ## build an x86_64 version of the docker container
	docker build --platform=linux/amd64 . -t cmaq

.phony: build-aarch64
build-aarch64: Dockerfile  ## build an arm version of the docker container
	docker build --platform=linux/arm64 . -t cmaq

.phony: build-conda
build-aarch64: Dockerfile  ## build just the "conda" step of the docker container
	docker build --platform=linux/arm64 . --target conda -t cmaq-conda

.phony: run
run: build  ## run the docker container
	docker run -it --rm -v ${PWD}:/opt/project cmaq


.phony: run-conda
run-conda: build  ## run a container with only the conda dependencies
	docker run -it --rm -v ${PWD}:/opt/project cmaq-conda