###################################################
#
# Makefile for gcat-lmm
#
###################################################

#
# Macros
#

# Header file locations for gcat-project
INCLUDE = -I/usr/include/gcat -I./src
# C++ compiler
CC = g++
# C++ compiler for MPICH2
MPICC = mpic++
# C++ linker
LD = g++
# C++ standard compiler options
CXXFLAGS = -Wall -w -O3 -g -D __NOEXTERN_FOR_CINCLUDE -D _NOT_MYUTILS_DEBUG
# C++ compiler options for gcat library code
CC_OPTIONS = $(CXXFLAGS) -fPIC
# C++ linker options
LNK_OPTIONS = -Wl,--no-as-needed -lxerces-c -L./

#
# Build gcat-project
#

LIB_GCAT_LMM_OBJECTS = \
		./Distribution_ContinuousVectorMixture.o\
		./Distribution_MarginalNormal.o\
		./Distribution_MultivariateNormal.o\
		./RandomVariable_ContinuousMatrix.o\
		./RandomVariable_ContinuousSymmetricMatrix.o\
		./Transformation_Centre.o\
		./Transformation_ComputeKinshipMatrix.o\
		./Transformation_ComputeKinshipMatrixPosition.o\
		./Transformation_ContinuousVectorMixtureComponentLogLikelihood.o\
		./Transformation_MeanVector.o\
		./Transformation_RescaleSymmetricMatrix.o\
		./Transformation_RescaleVector.o\
		./Transformation_SymmetricMatrixSum.o\
		./lmmLibrary.o\
		./lmmXML.o

all : libgcat_lmm.so

clean : cleanobj
	rm libgcat_lmm.so

cleanobj :
	rm $(LIB_GCAT_LMM_OBJECTS)

#install : gcat
#		cp gcat gcat

libgcat_lmm.so : $(LIB_GCAT_LMM_OBJECTS)
	$(LD) $(LNK_OPTIONS) -lgsl -lgcat-core -shared -o libgcat_lmm.so $(LIB_GCAT_LMM_OBJECTS)

#
# Build the parts of LMM
#

./Distribution_ContinuousVectorMixture.o : src/lmm/Distributions/ContinuousVectorMixture.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Distributions/ContinuousVectorMixture.cpp -c $(INCLUDE) -o ./Distribution_ContinuousVectorMixture.o

./Distribution_MarginalNormal.o : src/lmm/Distributions/MarginalNormal.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Distributions/MarginalNormal.cpp -c $(INCLUDE) -o ./Distribution_MarginalNormal.o

./Distribution_MultivariateNormal.o : src/lmm/Distributions/MultivariateNormal.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Distributions/MultivariateNormal.cpp -c $(INCLUDE) -o ./Distribution_MultivariateNormal.o

./RandomVariable_ContinuousMatrix.o : src/lmm/RandomVariables/ContinuousMatrix.cpp
	$(CC) $(CC_OPTIONS) src/lmm/RandomVariables/ContinuousMatrix.cpp -c $(INCLUDE) -o ./RandomVariable_ContinuousMatrix.o

./RandomVariable_ContinuousSymmetricMatrix.o : src/lmm/RandomVariables/ContinuousSymmetricMatrix.cpp
	$(CC) $(CC_OPTIONS) src/lmm/RandomVariables/ContinuousSymmetricMatrix.cpp -c $(INCLUDE) -o ./RandomVariable_ContinuousSymmetricMatrix.o

./Transformation_Centre.o : src/lmm/Transformations/Centre.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/Centre.cpp -c $(INCLUDE) -o ./Transformation_Centre.o

./Transformation_ComputeKinshipMatrix.o : src/lmm/Transformations/ComputeKinshipMatrix.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/ComputeKinshipMatrix.cpp -c $(INCLUDE) -o ./Transformation_ComputeKinshipMatrix.o

./Transformation_ComputeKinshipMatrixPosition.o : src/lmm/Transformations/ComputeKinshipMatrixPosition.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/ComputeKinshipMatrixPosition.cpp -c $(INCLUDE) -o ./Transformation_ComputeKinshipMatrixPosition.o

./Transformation_ContinuousVectorMixtureComponentLogLikelihood.o : src/lmm/Transformations/ContinuousVectorMixtureComponentLogLikelihood.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/ContinuousVectorMixtureComponentLogLikelihood.cpp -c $(INCLUDE) -o ./Transformation_ContinuousVectorMixtureComponentLogLikelihood.o

./Transformation_MeanVector.o : src/lmm/Transformations/MeanVector.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/MeanVector.cpp -c $(INCLUDE) -o ./Transformation_MeanVector.o

./Transformation_RescaleSymmetricMatrix.o : src/lmm/Transformations/RescaleSymmetricMatrix.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/RescaleSymmetricMatrix.cpp -c $(INCLUDE) -o ./Transformation_RescaleSymmetricMatrix.o

./Transformation_RescaleVector.o : src/lmm/Transformations/RescaleVector.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/RescaleVector.cpp -c $(INCLUDE) -o ./Transformation_RescaleVector.o

./Transformation_SymmetricMatrixSum.o : src/lmm/Transformations/SymmetricMatrixSumTransform.cpp
	$(CC) $(CC_OPTIONS) src/lmm/Transformations/SymmetricMatrixSumTransform.cpp -c $(INCLUDE) -o ./Transformation_SymmetricMatrixSum.o

./lmmLibrary.o : src/lmm/lmmLibrary.cpp
	$(CC) $(CC_OPTIONS) src/lmm/lmmLibrary.cpp -c $(INCLUDE) -o ./lmmLibrary.o

./lmmXML.o : src/lmm/lmmXML.cpp src/lmm/gcat_lmm1.0.xsd.h
	$(CC) $(CC_OPTIONS) src/lmm/lmmXML.cpp -c $(INCLUDE) -o ./lmmXML.o

src/lmm/gcat_lmm1.0.xsd.h : src/lmm/gcat_lmm1.0.xsd
	(cd src/lmm && xxd -i gcat_lmm1.0.xsd > gcat_lmm1.0.xsd.h)

##### END RUN ####
