all: $(addprefix bin/,bam_chunker_cascade mrsfast mrcanavar)

bin/bam_chunker_cascade:
	git submodule update --init --recursive dist/mrssnake
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
