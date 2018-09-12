$(shell mkdir -p bin)
OPENMP_CXX=g++ -lstdc++ -fopenmp -lz -lm

all: $(addprefix bin/,bam_chunker_cascade mrsfast mrcanavar samblaster)

bin/samblaster :
	git submodule update --init --recursive dist/samblaster
	-cd dist/samblaster && $(MAKE)
	-ln -s ../dist/samblaster/bin/samblaster bin/samblaster

bin/bam_chunker_cascade:
	git submodule update --init --recursive dist/mrssnake
	-cd dist/mrssnake && $(MAKE) CXX="$(OPENMP_CXX)"
	-@ln -s ../dist/mrssnake/bin/bam_chunker_cascade bin/bam_chunker_cascade

bin/mrsfast:
	git submodule update --init --recursive dist/mrsfast
	-cd dist/mrsfast && $(MAKE)
	-@ln -s ../dist/mrsfast/mrsfast bin/mrsfast
	-@ln -s ../dist/mrsfast/snp_indexer bin/snp_indexer

bin/mrcanavar:
	git submodule update --init --recursive dist/mrcanavar
	-cd dist/mrcanavar && $(MAKE)
	-@ln -s ../dist/mrcanavar/mrcanavar bin/mrcanavar
