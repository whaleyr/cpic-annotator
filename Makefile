.ONESHELL:
.SHELLFLAGS = -e -c
MAKEFLAGS = -j1

# Determine this makefile's path.
# Be sure to place this BEFORE `include` directives, if any.
THIS_FILE := $(lastword $(MAKEFILE_LIST))
ROOT_DIR := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

ifdef PHARMCAT_DATA_DIR
  dataDir := $(PHARMCAT_DATA_DIR)
else
  dataDir := build
endif

ifeq ($(OS),Windows_NT)
	GRADLE_CMD := cmd /c gradlew.bat --console=plain
	ifneq (,$(findstring cygdrive,${ROOT_DIR}))
		ROOT_DIR := $(shell echo ${ROOT_DIR} | sed 's/\/cygdrive\///' | sed 's/./&:/1')
	endif
else
	GRADLE_CMD := ./gradlew --console=plain
endif

ifdef PHARMCAT_DOCKER_PLATFORM
	dockerPlatform := --platform $(PHARMCAT_DOCKER_PLATFORM)
else
	dockerPlatform :=
endif

quietTest = true


.PHONY: updateData
updateData: clean
	@echo "Updating data..."
	@${GRADLE_CMD} updateData
	@echo "Moving pharmcat_position files..."
	mv -f src/main/resources/org/pharmgkb/pharmcat/definition/alleles/pharmcat_positions.* .
	cp -f pharmcat_positions.vcf src/test/resources/org/pharmgkb/pharmcat/reference.vcf
	@echo ""
	@echo "Updating examples..."
	@${GRADLE_CMD} updateExample
	# this reverts files with only EOL changes
	@git stash
	@git stash pop


.PHONY: _dockerPrep
_dockerPrep: clean
	${GRADLE_CMD} shadowJar
	mv -f build/libs/pharmcat-`git describe --tags | sed -r s/^v//`-all.jar build/pharmcat.jar

.PHONY: docker
docker: _dockerPrep
	docker build ${dockerPlatform} -t pcat .


.PHONY: scriptPkg
scriptPkg:
	rm -rf build/preprocessor
	mkdir -p build/preprocessor
	cp -f preprocessor/requirements.txt build/preprocessor
	cp -f preprocessor/pharmcat_vcf_preprocessor.py build/preprocessor
	cp -f preprocessor/pharmcat_pipeline build/preprocessor
	mkdir build/preprocessor/preprocessor
	cp -f preprocessor/preprocessor/*.py build/preprocessor/preprocessor
	cp -f preprocessor/preprocessor/*.tsv build/preprocessor/preprocessor
	cp -f pharmcat_positions.vcf* build/preprocessor
	cp -f docs/using/VCF-Preprocessor.md build/preprocessor/README.md
	cd build; tar -czvf preprocessor.tar.gz preprocessor


.PHONE: updateDataFromScratch
updateDataFromScratch: docker updateData


.PHONY: clean
clean:
	@echo "Cleaning up..."
	@${GRADLE_CMD} -q clean

.PHONY: clean-test-data
clean-test-data:
	@echo "Deleting pre-existing VCF test data"
	@rm -rf ${dataDir}/testVcf

.PHONY: clean-test-zips
clean-test-zips:
	@echo "Deleting pre-existing VCF test zips"
	rm -f ${dataDir}/vcfTests-*.zip


.PHONY: allVcfTests
allVcfTests: clean _generateVcf _fuzzyVcfTests _vcfTests _exactVcfTests

.PHONY: rerun-allVcfTests
rerun-allVcfTests: clean _fuzzyVcfTests _vcfTests _exactVcfTests

.PHONY: allVcfMissingTests
allVcfMissingTests: clean _generateVcf-missing _fuzzyVcfTests-missing _vcfTests-missing _exactVcfTests-missing

.PHONY: _generateVcf
_generateVcf: clean-test-data
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Generating VCF test data"
	@src/scripts/vcf_generator/generate_vcf_test_data.sh
	@echo ""

.PHONY: _generateVcf-missing
_generateVcf-missing: clean-test-data
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Generating VCF test data with missing positions"
	@src/scripts/vcf_generator/generate_vcf_test_data.sh -m
	@echo ""


.PHONY: _vcfTests
_vcfTests:
	@echo "Deleting pre-existing VCF test results"
	@rm -rf ${dataDir}/autogeneratedTestResults
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Running default VCF tests:"
	@${GRADLE_CMD} -q testAutogeneratedVcfs
	@echo "Zipping results..."
	@src/scripts/vcf_generator/zip_results.sh
	@src/scripts/vcf_generator/compare_results.sh

.PHONY: vcfTests
vcfTests: clean _generateVcf _vcfTests


.PHONY: _vcfTests-missing
_vcfTests-missing:
	@echo "Deleting pre-existing VCF test results"
	@rm -rf ${dataDir}/autogeneratedTestResults
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Running default VCF with missing positions tests:"
	@${GRADLE_CMD} -q testAutogeneratedVcfsMissing
	@echo "Zipping results..."
	@src/scripts/vcf_generator/zip_results.sh -m
	@src/scripts/vcf_generator/compare_results.sh -m

.PHONY: vcfMissingTests
vcfMissingTests: clean _generateVcf-missing _vcfTests-missing


.PHONY: _exactVcfTests
_exactVcfTests:
	@echo "Deleting pre-existing VCF test results"
	@rm -rf ${dataDir}/autogeneratedTestResults
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Running exact VCF tests:"
	@${GRADLE_CMD} -q testAutogeneratedVcfsExactMatchOnly
	@echo "Zipping results..."
	src/scripts/vcf_generator/zip_results.sh -e
	src/scripts/vcf_generator/compare_results.sh -e

.PHONY: exactVcfTests
exactVcfTests: clean _generateVcf _exactVcfTests


.PHONY: _exactVcfTests-missing
_exactVcfTests-missing:
	@echo "Deleting pre-existing VCF test results"
	@rm -rf ${dataDir}/autogeneratedTestResults
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Running exact VCF with missing positions tests:"
	@${GRADLE_CMD} -q testAutogeneratedVcfsExactMatchOnlyMissing
	@echo "Zipping results..."
	src/scripts/vcf_generator/zip_results.sh -e -m
	src/scripts/vcf_generator/compare_results.sh -e -m

.PHONY: exactVcfMissingTests
exactVcfMissingTests: clean _generateVcf-missing _exactVcfTests-missing


.PHONY: _fuzzyVcfTests
_fuzzyVcfTests:
	@echo "Deleting pre-existing VCF test results"
	@rm -rf ${dataDir}/autogeneratedTestResults
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Running fuzzy VCF tests:"
	${GRADLE_CMD} -q testAutogeneratedVcfsFuzzyMatch
	@echo "Zipping results..."
	src/scripts/vcf_generator/zip_results.sh -f
	src/scripts/vcf_generator/compare_results.sh -f

.PHONY: fuzzyVcfTests
fuzzyVcfTests: clean _generateVcf _fuzzyVcfTests


.PHONY: _fuzzyVcfTests-missing
_fuzzyVcfTests-missing:
	@echo "Deleting pre-existing VCF test results"
	@rm -rf ${dataDir}/autogeneratedTestResults
	@export PHARMCAT_TEST_QUIET=${quietTest}
	@echo "Running fuzzy VCF with missing positions tests:"
	${GRADLE_CMD} -q testAutogeneratedVcfsFuzzyMatchMissing
	@echo "Zipping results..."
	src/scripts/vcf_generator/zip_results.sh -f -m
	src/scripts/vcf_generator/compare_results.sh -f -m

.PHONY: fuzzyVcfMissingTests
fuzzyVcfMissingTests: clean _generateVcf-missing _fuzzyVcfTests-missing


.PHONY: test-preprocessor
test-preprocessor:
	cd preprocessor
	python3 -m pytest --cov=preprocessor


.PHONY: release
release:
	bin/release.sh


.PHONY: dockerRelease
dockerRelease: _dockerPrep
	version=`git describe --tags | sed -r s/^v//`
	docker buildx inspect pcat-builder || docker buildx create --name pcat-builder --bootstrap --platform linux/amd64,linux/arm64
	docker buildx build --builder pcat-builder --platform linux/amd64,linux/arm64 --push --tag pgkb/pharmcat:latest --tag pgkb/pharmcat:$${version} .

.PHONEY: jekyllDocker
jekyllDocker:
	docker build ${dockerPlatform} -t jekyll docs/

.PHONEY: jekyll
jekyll:
	docker run --rm -it -v ${ROOT_DIR}:/PharmCAT jekyll


# used to build reference FASTA download for Zenodo
.PHONEY: referenceFasta
referenceFasta:
	bin/referenceFasta.sh
