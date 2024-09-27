TAG = weisburd/tandem-repeat-catalogs:latest

all: build push

build:
	cd .. && docker build -t $(TAG) -f docker/Dockerfile .

push:
	docker push $(TAG)
	docker push $(TAG) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
